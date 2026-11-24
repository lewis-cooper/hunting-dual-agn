"""Pair geometry, Gaia association and fixed-primary assignment from the existing search."""

import math
from pathlib import Path
from typing import Any, Mapping

from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table, join
import numpy as np
import pandas as pd


def finite_float(value):
    if np.ma.is_masked(value):
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def text_or_empty(value):
    return '' if pd.isna(value) else str(value).strip()


def require_gaia_coverage(components, coverage, radius_arcsec):
    required = {'ra_deg', 'dec_deg', 'query_radius_arcsec'}
    if not required.issubset(coverage):
        raise ValueError('Gaia coverage needs ra_deg, dec_deg, query_radius_arcsec for every successful component cone, including empty results.')
    valid = np.isfinite(coverage[['ra_deg','dec_deg','query_radius_arcsec']].to_numpy(float)).all(axis=1)
    coverage = coverage.loc[valid & coverage.query_radius_arcsec.ge(radius_arcsec)].copy()
    if coverage.empty:
        raise ValueError('No valid Gaia component-cone coverage was supplied')
    positions = SkyCoord(np.asarray(components['component_ra'])*u.deg,
                         np.asarray(components['component_dec'])*u.deg)
    contained = np.zeros(len(components),dtype=bool)
    for covered_radius, group in coverage.groupby('query_radius_arcsec'):
        centres = SkyCoord(group.ra_deg.to_numpy()*u.deg, group.dec_deg.to_numpy()*u.deg)
        _, offsets, _ = positions.match_to_catalog_sky(centres)
        contained |= offsets.arcsec + radius_arcsec <= covered_radius + 1e-6
    if not contained.all():
        raise ValueError(f'{(~contained).sum()} component cones lack complete Gaia coverage. Supply local data; no download was attempted.')


def match_gaia_components(components, gaia, radius_arcsec):
    """Match against local Gaia sources; non-detections require separate coverage proof."""
    required = ['source_id','ra','dec','phot_g_mean_mag','pmra','pmdec','pmra_error',
                'pmdec_error','pmra_pmdec_corr','parallax_over_error','ruwe']
    missing = set(required).difference(gaia.columns)
    if missing:
        raise ValueError('Gaia catalogue missing '+', '.join(sorted(missing)))
    gaia = gaia[required].drop_duplicates('source_id').copy()
    if not gaia.source_id.astype('string').str.fullmatch('[1-9][0-9]*').all():
        raise ValueError('Gaia source_id must be an exact positive identifier')
    gaia['source_id'] = gaia.source_id.astype('int64')
    gaia['pm_sig'] = pm_sig(gaia.pmra,gaia.pmdec,gaia.pmra_error,gaia.pmdec_error,gaia.pmra_pmdec_corr)
    component_sky = SkyCoord(np.asarray(components['component_ra'])*u.deg,
                             np.asarray(components['component_dec'])*u.deg)
    gaia_sky = SkyCoord(gaia.ra.to_numpy()*u.deg,gaia.dec.to_numpy()*u.deg)
    component_index, gaia_index, offsets, _ = search_around_sky(component_sky,gaia_sky,radius_arcsec*u.arcsec)
    matches = gaia.iloc[gaia_index].reset_index(drop=True)
    matches['component_row_index'] = component_index
    matches['component_gaia_sep_arcsec'] = offsets.arcsec
    matches = matches.sort_values(['component_row_index','component_gaia_sep_arcsec'],kind='stable')
    return Table.from_pandas(matches.reset_index(drop=True))



def angular_separation_arcsec(
    ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float
) -> float:
    delta_ra = ((ra1_deg - ra2_deg + 180.0) % 360.0) - 180.0
    mean_dec_rad = math.radians(0.5 * (dec1_deg + dec2_deg))
    return 3600.0 * math.hypot(delta_ra * math.cos(mean_dec_rad), dec1_deg - dec2_deg)


def flux_colour(numerator_flux: object, denominator_flux: object) -> float:
    numerator = finite_float(numerator_flux)
    denominator = finite_float(denominator_flux)
    if not np.isfinite(numerator) or not np.isfinite(denominator):
        return math.nan
    if numerator <= 0 or denominator <= 0:
        return math.nan
    return float(-2.5 * np.log10(numerator / denominator))


def magnitude_from_flux(flux_value: object) -> float:
    flux = finite_float(flux_value)
    if not np.isfinite(flux) or flux <= 0:
        return math.nan
    return float(22.5 - 2.5 * np.log10(flux))


def channel_row(
    parent: Mapping[str, object], local: pd.DataFrame, aperture_arcsec: float
) -> tuple[dict[str, object] | None, pd.DataFrame]:
    if local.empty:
        return None, local.iloc[0:0]
    separation = pd.to_numeric(local["target_sep_arcsec"], errors="coerce").to_numpy(float)
    flux_r = pd.to_numeric(local["flux_r"], errors="coerce").to_numpy(float)
    in_aperture = np.flatnonzero(
        np.isfinite(separation) & (separation <= aperture_arcsec) & np.isfinite(flux_r)
    )
    if in_aperture.size < 2:
        return None, local.iloc[0:0]
    ranked = in_aperture[np.argsort(flux_r[in_aperture], kind="stable")[::-1]]
    primary = local.iloc[int(ranked[0])]
    secondary = local.iloc[int(ranked[1])]
    primary_flux_r = finite_float(primary.get("flux_r"))
    secondary_flux_r = finite_float(secondary.get("flux_r"))
    primary_g_minus_r = flux_colour(primary.get("flux_g"), primary_flux_r)
    secondary_g_minus_r = flux_colour(secondary.get("flux_g"), secondary_flux_r)
    primary_r_minus_z = flux_colour(primary_flux_r, primary.get("flux_z"))
    secondary_r_minus_z = flux_colour(secondary_flux_r, secondary.get("flux_z"))
    primary_type = text_or_empty(primary.get("type"))
    secondary_type = text_or_empty(secondary.get("type"))
    if primary_type == "PSF" and secondary_type == "PSF":
        type_pair = "PSF-PSF"
    elif "PSF" in {primary_type, secondary_type}:
        type_pair = "PSF-extended"
    else:
        type_pair = "non-PSF-pair"
    row = dict(parent)
    row.update(
        {
            "n_within_aperture": int(in_aperture.size),
            "primary_type": primary_type,
            "secondary_type": secondary_type,
            "type_pair": type_pair,
            "primary_target_sep_arcsec": finite_float(primary.get("target_sep_arcsec")),
            "secondary_target_sep_arcsec": finite_float(secondary.get("target_sep_arcsec")),
            "primary_ra": finite_float(primary.get("ra")),
            "primary_dec": finite_float(primary.get("dec")),
            "secondary_ra": finite_float(secondary.get("ra")),
            "secondary_dec": finite_float(secondary.get("dec")),
            "component_sep_arcsec": angular_separation_arcsec(
                finite_float(primary.get("ra")),
                finite_float(primary.get("dec")),
                finite_float(secondary.get("ra")),
                finite_float(secondary.get("dec")),
            ),
            "primary_flux_r": primary_flux_r,
            "secondary_flux_r": secondary_flux_r,
            "secondary_to_primary_flux_r": secondary_flux_r / primary_flux_r
            if np.isfinite(primary_flux_r)
            and primary_flux_r > 0
            and np.isfinite(secondary_flux_r)
            else math.nan,
            "primary_r_mag": magnitude_from_flux(primary_flux_r),
            "secondary_r_mag": magnitude_from_flux(secondary_flux_r),
            "primary_g_minus_r": primary_g_minus_r,
            "secondary_g_minus_r": secondary_g_minus_r,
            "delta_g_minus_r": secondary_g_minus_r - primary_g_minus_r
            if np.isfinite(primary_g_minus_r) and np.isfinite(secondary_g_minus_r)
            else math.nan,
            "primary_r_minus_z": primary_r_minus_z,
            "secondary_r_minus_z": secondary_r_minus_z,
            "delta_r_minus_z": secondary_r_minus_z - primary_r_minus_z
            if np.isfinite(primary_r_minus_z) and np.isfinite(secondary_r_minus_z)
            else math.nan,
        }
    )
    return row, local.copy()


def float_array(table: Table, name: str) -> np.ndarray:
    return np.asarray(table[name], dtype=float)


def str_array(table: Table, name: str) -> np.ndarray:
    return np.char.strip(np.asarray(table[name]).astype(str))


def integer_scalar_or_default(value: object, default: int = -1) -> int:
    """Parse an integer scalar without losing precision for string Gaia IDs."""
    if value is None or np.ma.is_masked(value):
        return default
    if isinstance(value, (bytes, np.bytes_)):
        value = value.decode("ascii", errors="strict")
    if isinstance(value, (str, np.str_)):
        text = str(value).strip()
        if not text:
            return default
        try:
            return int(text)
        except ValueError:
            value = float(text)
    try:
        numeric_value = float(value)
    except (TypeError, ValueError, OverflowError):
        return default
    if not np.isfinite(numeric_value):
        return default
    return int(value)


def component_sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    dra = ((ra1 - ra2 + 180.0) % 360.0) - 180.0
    dec_mean = np.deg2rad(0.5 * (dec1 + dec2))
    return float(3600.0 * np.hypot(dra * np.cos(dec_mean), dec1 - dec2))


def component_roles(components_per_system: int) -> tuple[str, ...]:
    if not 2 <= int(components_per_system) <= 6:
        raise ValueError("components_per_system must be between 2 and 6")
    base_roles = ("primary", "secondary", "tertiary")
    extra_roles = tuple(f"component{rank}" for rank in range(4, int(components_per_system) + 1))
    return base_roles[: min(components_per_system, len(base_roles))] + extra_roles


def build_component_table(
    targets: Table,
    tractor: Table,
    aperture_arcsec: float,
    max_targets: int,
    components_per_system: int = 2,
    min_components_per_system: int = 2,
) -> Table:
    roles = component_roles(components_per_system)
    primary_type = str_array(targets, "primary_type")
    secondary_type = str_array(targets, "secondary_type")
    ch1a = targets[(primary_type == "PSF") & (secondary_type == "PSF")]
    if max_targets > 0:
        ch1a = ch1a[: int(max_targets)]

    selected_qidx = set(int(value) for value in np.asarray(ch1a["qso_parent_index"], dtype=np.int64))
    tractor_qidx = np.asarray(tractor["qso_parent_index"], dtype=np.int64)
    keep = np.array([int(value) in selected_qidx for value in tractor_qidx], dtype=bool)
    sub = tractor[keep]

    order = np.argsort(np.asarray(sub["qso_parent_index"], dtype=np.int64), kind="stable")
    qidx_sorted = np.asarray(sub["qso_parent_index"], dtype=np.int64)[order]
    unique_qidx, starts, counts = np.unique(qidx_sorted, return_index=True, return_counts=True)

    target_by_qidx = {int(row["qso_parent_index"]): row for row in ch1a}
    sep = float_array(sub, "target_sep_arcsec")[order]
    typ = str_array(sub, "type")[order]
    flux_r = float_array(sub, "flux_r")[order]
    flux_g = float_array(sub, "flux_g")[order]
    flux_z = float_array(sub, "flux_z")[order]
    ra = float_array(sub, "ra")[order]
    dec = float_array(sub, "dec")[order]
    objid = np.asarray(sub["objid"])[order]
    brickid = np.asarray(sub["brickid"])[order]
    release = np.asarray(sub["release"])[order]

    rows: list[dict[str, object]] = []
    for qidx, start, count in zip(unique_qidx, starts, counts):
        qidx_int = int(qidx)
        target = target_by_qidx.get(qidx_int)
        if target is None:
            continue
        sl = slice(int(start), int(start + count))
        in_ap = np.flatnonzero((sep[sl] <= aperture_arcsec) & (typ[sl] == "PSF"))
        if len(in_ap) < min_components_per_system:
            continue
        group = int(start) + in_ap
        rank_flux = flux_r[group]
        finite_rank = np.isfinite(rank_flux) & (rank_flux > 0)
        if np.sum(finite_rank) < min_components_per_system:
            continue
        finite_group = group[finite_rank]
        ranked = finite_group[np.argsort(flux_r[finite_group])[::-1]][:components_per_system]
        comp_indices = {role: int(idx) for role, idx in zip(roles, ranked)}
        primary = comp_indices["primary"]
        secondary = comp_indices["secondary"]
        tertiary = comp_indices.get("tertiary")
        comp_sep = component_sep_arcsec(ra[primary], dec[primary], ra[secondary], dec[secondary])
        primary_tertiary_sep = (
            component_sep_arcsec(ra[primary], dec[primary], ra[tertiary], dec[tertiary]) if tertiary is not None else math.nan
        )
        secondary_tertiary_sep = (
            component_sep_arcsec(ra[secondary], dec[secondary], ra[tertiary], dec[tertiary]) if tertiary is not None else math.nan
        )
        g_minus_r = {role: flux_colour(flux_g[idx], flux_r[idx]) for role, idx in comp_indices.items()}
        r_minus_z = {role: flux_colour(flux_r[idx], flux_z[idx]) for role, idx in comp_indices.items()}
        p_gr = g_minus_r["primary"]
        s_gr = g_minus_r["secondary"]
        p_rz = r_minus_z["primary"]
        s_rz = r_minus_z["secondary"]
        colour_values = np.asarray(list(g_minus_r.values()), dtype=float)
        finite_colours = colour_values[np.isfinite(colour_values)]
        max_pairwise_delta_gr = (
            float(np.nanmax(finite_colours) - np.nanmin(finite_colours)) if finite_colours.size == len(comp_indices) else math.nan
        )
        ratio = (
            float(flux_r[secondary] / flux_r[primary])
            if np.isfinite(flux_r[primary]) and flux_r[primary] > 0 and np.isfinite(flux_r[secondary])
            else math.nan
        )
        for role, idx in comp_indices.items():
            component_delta_gr = (
                float(g_minus_r[role] - p_gr) if np.isfinite(g_minus_r[role]) and np.isfinite(p_gr) else math.nan
            )
            component_delta_rz = (
                float(r_minus_z[role] - p_rz) if np.isfinite(r_minus_z[role]) and np.isfinite(p_rz) else math.nan
            )
            rows.append(
                {
                    "component_id": f"{qidx_int}_{role}",
                    "qso_parent_index": qidx_int,
                    "sdss_name": str(target["sdss_name"]),
                    "z": float(target["z"]),
                    "component_role": role,
                    "component_rank_r": int(roles.index(role) + 1),
                    "n_psf_components_within_aperture": int(len(in_ap)),
                    "component_ra": float(ra[idx]),
                    "component_dec": float(dec[idx]),
                    "component_target_sep_arcsec": float(sep[idx]),
                    "component_sep_arcsec": comp_sep,
                    "primary_secondary_sep_arcsec": comp_sep,
                    "primary_tertiary_sep_arcsec": primary_tertiary_sep,
                    "secondary_tertiary_sep_arcsec": secondary_tertiary_sep,
                    "tractor_release": int(release[idx]),
                    "tractor_brickid": int(brickid[idx]),
                    "tractor_objid": int(objid[idx]),
                    "tractor_flux_g": float(flux_g[idx]),
                    "tractor_flux_r": float(flux_r[idx]),
                    "tractor_flux_z": float(flux_z[idx]),
                    "tractor_r_mag": magnitude_from_flux(float(flux_r[idx])),
                    "tractor_g_minus_r": g_minus_r[role],
                    "tractor_r_minus_z": r_minus_z[role],
                    "secondary_to_primary_flux_r": ratio,
                    "delta_r_mag": -2.5 * np.log10(ratio) if np.isfinite(ratio) and ratio > 0 else math.nan,
                    "delta_g_minus_r": float(s_gr - p_gr) if np.isfinite(p_gr) and np.isfinite(s_gr) else math.nan,
                    "delta_r_minus_z": float(s_rz - p_rz) if np.isfinite(p_rz) and np.isfinite(s_rz) else math.nan,
                    "component_delta_g_minus_r_from_primary": component_delta_gr,
                    "component_delta_r_minus_z_from_primary": component_delta_rz,
                    "max_pairwise_abs_delta_g_minus_r": max_pairwise_delta_gr,
                    "original_gaia_source_id": integer_scalar_or_default(target["source_id"]),
                    "original_gaia_ra": float(target["gaia_ra"]),
                    "original_gaia_dec": float(target["gaia_dec"]),
                    "q_ra": float(target["q_ra"]),
                    "q_dec": float(target["q_dec"]),
                }
            )
    out = Table(rows=rows)
    out.sort(["qso_parent_index", "component_rank_r"])
    return out


def pm_sig(pmra: np.ndarray, pmdec: np.ndarray, pmra_error: np.ndarray, pmdec_error: np.ndarray, corr: np.ndarray) -> np.ndarray:
    pmra = np.asarray(pmra, dtype=float)
    pmdec = np.asarray(pmdec, dtype=float)
    era = np.asarray(pmra_error, dtype=float)
    edec = np.asarray(pmdec_error, dtype=float)
    rho = np.asarray(corr, dtype=float)
    out = np.full(pmra.shape, np.nan, dtype=float)
    valid = (
        np.isfinite(pmra)
        & np.isfinite(pmdec)
        & np.isfinite(era)
        & np.isfinite(edec)
        & np.isfinite(rho)
        & (era > 0)
        & (edec > 0)
        & (np.abs(rho) < 1)
    )
    chi2 = (
        (pmra[valid] / era[valid]) ** 2
        + (pmdec[valid] / edec[valid]) ** 2
        - 2.0 * rho[valid] * pmra[valid] * pmdec[valid] / (era[valid] * edec[valid])
    ) / (1.0 - rho[valid] ** 2)
    out[valid] = np.sqrt(np.maximum(chi2, 0.0))
    return out


def first_match_by_component(matches: Table) -> dict[int, Any]:
    out: dict[int, Any] = {}
    if not len(matches):
        return out
    for i, row in enumerate(matches):
        comp_idx = int(row["component_row_index"])
        if comp_idx not in out:
            out[comp_idx] = row
    return out


def build_candidate_summary(components: Table, matches: Table) -> Table:
    first = first_match_by_component(matches)
    component_rows = {int(i): row for i, row in enumerate(components)}
    rows: list[dict[str, object]] = []
    for qidx in np.unique(np.asarray(components["qso_parent_index"], dtype=np.int64)):
        qidx_int = int(qidx)
        comp_idx = np.flatnonzero(np.asarray(components["qso_parent_index"], dtype=np.int64) == qidx_int)
        if len(comp_idx) < 2:
            continue
        comp_idx = comp_idx[np.argsort(np.asarray(components["component_rank_r"][comp_idx], dtype=int))]
        p_i, s_i = int(comp_idx[0]), int(comp_idx[1])
        p = component_rows[p_i]
        s = component_rows[s_i]
        t_i = int(comp_idx[2]) if len(comp_idx) >= 3 else None
        t = component_rows[t_i] if t_i is not None else None
        p_match = first.get(p_i)
        s_match = first.get(s_i)
        t_match = first.get(t_i) if t_i is not None else None

        p_sid = int(p_match["source_id"]) if p_match is not None else -1
        s_sid = int(s_match["source_id"]) if s_match is not None else -1
        t_sid = int(t_match["source_id"]) if t_match is not None else -1
        original_sid = int(p["original_gaia_source_id"])
        row = {
            "qso_parent_index": qidx_int,
            "sdss_name": str(p["sdss_name"]),
            "z": finite_float(p["z"]),
            "component_sep_arcsec": finite_float(p["component_sep_arcsec"]),
            "primary_secondary_sep_arcsec": finite_float(p["primary_secondary_sep_arcsec"]),
            "primary_tertiary_sep_arcsec": float(np.ma.filled(p["primary_tertiary_sep_arcsec"],math.nan)) if "primary_tertiary_sep_arcsec" in p.colnames else math.nan,
            "secondary_tertiary_sep_arcsec": float(np.ma.filled(p["secondary_tertiary_sep_arcsec"],math.nan)) if "secondary_tertiary_sep_arcsec" in p.colnames else math.nan,
            "n_psf_components_within_aperture": int(p["n_psf_components_within_aperture"])
            if "n_psf_components_within_aperture" in p.colnames
            else len(comp_idx),
            "secondary_to_primary_flux_r": finite_float(p["secondary_to_primary_flux_r"]),
            "delta_r_mag": finite_float(p["delta_r_mag"]),
            "delta_g_minus_r": finite_float(p["delta_g_minus_r"]),
            "delta_r_minus_z": finite_float(p["delta_r_minus_z"]) if "delta_r_minus_z" in p.colnames else math.nan,
            "max_pairwise_abs_delta_g_minus_r": finite_float(p["max_pairwise_abs_delta_g_minus_r"])
            if "max_pairwise_abs_delta_g_minus_r" in p.colnames
            else math.nan,
            "primary_ra": finite_float(p["component_ra"]),
            "primary_dec": finite_float(p["component_dec"]),
            "secondary_ra": finite_float(s["component_ra"]),
            "secondary_dec": finite_float(s["component_dec"]),
            "original_gaia_source_id": original_sid,
            "primary_gaia_source_id": p_sid,
            "primary_gaia_sep_arcsec": finite_float(p_match["component_gaia_sep_arcsec"]) if p_match is not None else math.nan,
            "primary_gaia_g": finite_float(p_match["phot_g_mean_mag"]) if p_match is not None else math.nan,
            "primary_pm_sig": finite_float(p_match["pm_sig"]) if p_match is not None else math.nan,
            "primary_parallax_over_error": finite_float(p_match["parallax_over_error"]) if p_match is not None else math.nan,
            "primary_ruwe": finite_float(p_match["ruwe"]) if p_match is not None else math.nan,
            "secondary_gaia_source_id": s_sid,
            "secondary_gaia_sep_arcsec": finite_float(s_match["component_gaia_sep_arcsec"]) if s_match is not None else math.nan,
            "secondary_gaia_g": finite_float(s_match["phot_g_mean_mag"]) if s_match is not None else math.nan,
            "secondary_pm_sig": finite_float(s_match["pm_sig"]) if s_match is not None else math.nan,
            "secondary_parallax_over_error": finite_float(s_match["parallax_over_error"]) if s_match is not None else math.nan,
            "secondary_ruwe": finite_float(s_match["ruwe"]) if s_match is not None else math.nan,
            "primary_has_gaia": p_sid >= 0,
            "secondary_has_gaia": s_sid >= 0,
            "secondary_has_distinct_gaia": (s_sid >= 0) and (s_sid != p_sid),
            "secondary_is_original_gaia": (s_sid >= 0) and (s_sid == original_sid),
            "primary_is_original_gaia": (p_sid >= 0) and (p_sid == original_sid),
        }
        if t is not None:
            row.update(
                {
                    "tertiary_ra": float(t["component_ra"]),
                    "tertiary_dec": float(t["component_dec"]),
                    "tertiary_g_minus_r": float(t["tractor_g_minus_r"]),
                    "tertiary_r_minus_z": float(t["tractor_r_minus_z"]),
                    "tertiary_delta_g_minus_r": float(t["component_delta_g_minus_r_from_primary"]),
                    "tertiary_delta_r_minus_z": float(t["component_delta_r_minus_z_from_primary"]),
                    "tertiary_gaia_source_id": t_sid,
                    "tertiary_gaia_sep_arcsec": float(t_match["component_gaia_sep_arcsec"]) if t_match is not None else math.nan,
                    "tertiary_gaia_g": float(t_match["phot_g_mean_mag"]) if t_match is not None else math.nan,
                    "tertiary_pm_sig": float(t_match["pm_sig"]) if t_match is not None else math.nan,
                    "tertiary_parallax_over_error": float(t_match["parallax_over_error"]) if t_match is not None else math.nan,
                    "tertiary_ruwe": float(t_match["ruwe"]) if t_match is not None else math.nan,
                    "tertiary_has_gaia": t_sid >= 0,
                    "tertiary_has_distinct_gaia": (t_sid >= 0) and (t_sid not in {p_sid, s_sid}),
                    "tertiary_is_original_gaia": (t_sid >= 0) and (t_sid == original_sid),
                    "n_distinct_gaia_sources": len({sid for sid in (p_sid, s_sid, t_sid) if sid >= 0}),
                    "all_three_components_have_distinct_gaia": (p_sid >= 0) and (s_sid >= 0) and (t_sid >= 0) and len({p_sid, s_sid, t_sid}) == 3,
                }
            )
        rows.append(row)
    return Table(rows=rows)


def local_color_stats(
    query_qidx: np.ndarray,
    query_z: np.ndarray,
    ref_qidx: np.ndarray,
    ref_z: np.ndarray,
    ref_color: np.ndarray,
    z_tol: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(ref_z)
    z_sorted = ref_z[order]
    color_sorted = ref_color[order]
    qidx_sorted = ref_qidx[order]
    med = np.full(len(query_z), np.nan, dtype=float)
    sig = np.full(len(query_z), np.nan, dtype=float)
    p16 = np.full(len(query_z), np.nan, dtype=float)
    p84 = np.full(len(query_z), np.nan, dtype=float)
    n_ref = np.zeros(len(query_z), dtype=int)
    for i, (qidx, z) in enumerate(zip(query_qidx, query_z)):
        if not np.isfinite(z):
            continue
        lo = np.searchsorted(z_sorted, z - z_tol, side="left")
        hi = np.searchsorted(z_sorted, z + z_tol, side="right")
        vals = color_sorted[lo:hi]
        ids = qidx_sorted[lo:hi]
        vals = vals[ids != qidx]
        vals = vals[np.isfinite(vals)]
        if vals.size < 20:
            continue
        p16_i, med_i, p84_i = np.nanpercentile(vals, [16.0, 50.0, 84.0])
        spread = 0.5 * (p84_i - p16_i)
        if not np.isfinite(spread) or spread <= 0:
            continue
        p16[i] = p16_i
        med[i] = med_i
        p84[i] = p84_i
        sig[i] = spread
        n_ref[i] = vals.size
    return med, sig, p16, p84, n_ref


def as_float(table: Table, name: str) -> np.ndarray:
    return np.asarray(table[name], dtype=float)


def finite_source_id(table: Table, name: str) -> np.ndarray:
    if name not in table.colnames:
        return np.full(len(table), -1, dtype=np.int64)
    values = np.asarray(table[name], dtype=np.int64)
    return values


def swap_columns(table: Table, first: str, second: str, swap_mask: np.ndarray) -> None:
    if first not in table.colnames or second not in table.colnames:
        return
    first_values = np.array(table[first], copy=True)
    second_values = np.array(table[second], copy=True)
    table[first][swap_mask] = second_values[swap_mask]
    table[second][swap_mask] = first_values[swap_mask]


def relabel_by_qso_colour(
    summary: Table,
    targets: Table,
    *,
    reference: pd.DataFrame,
    z_tol: float,
) -> tuple[Table, dict[str, object]]:
    keep_cols = ["qso_parent_index", "primary_g_minus_r", "secondary_g_minus_r", "q_ra", "q_dec"]
    missing = sorted(set(keep_cols) - set(targets.colnames))
    if missing:
        raise ValueError(f"Target table is missing required columns: {missing}")

    joined = join(summary, targets[keep_cols], keys="qso_parent_index", join_type="left")
    qidx = np.asarray(joined["qso_parent_index"], dtype=np.int64)
    redshift = as_float(joined, "z")
    primary_gr = as_float(joined, "primary_g_minus_r")
    secondary_gr = as_float(joined, "secondary_g_minus_r")

    ref_qidx = reference.qso_parent_index.to_numpy()
    ref_z = reference.z.to_numpy()
    ref_color = reference.nearest_g_minus_r.to_numpy()
    med, sig, p16, p84, n_ref = local_color_stats(qidx, redshift, ref_qidx, ref_z, ref_color, z_tol)

    primary_scaled = np.full(len(joined), np.nan, dtype=float)
    secondary_scaled = np.full(len(joined), np.nan, dtype=float)
    finite_reference = np.isfinite(med) & np.isfinite(sig) & (sig > 0.0)
    primary_scaled[finite_reference] = np.abs((primary_gr[finite_reference] - med[finite_reference]) / sig[finite_reference])
    secondary_scaled[finite_reference] = np.abs((secondary_gr[finite_reference] - med[finite_reference]) / sig[finite_reference])

    primary_colour_ok = np.isfinite(primary_scaled)
    secondary_colour_ok = np.isfinite(secondary_scaled)
    can_use_colour = primary_colour_ok | secondary_colour_ok
    old_primary_is_qso_colour = (
        primary_colour_ok
        & (
            ~secondary_colour_ok
            | (primary_scaled <= secondary_scaled)
        )
    )

    old_primary_is_nearest = np.full(len(joined), True, dtype=bool)
    nearest_ok = np.zeros(len(joined), dtype=bool)
    coord_mask = (
        np.isfinite(as_float(joined, "q_ra"))
        & np.isfinite(as_float(joined, "q_dec"))
        & np.isfinite(as_float(joined, "primary_ra"))
        & np.isfinite(as_float(joined, "primary_dec"))
        & np.isfinite(as_float(joined, "secondary_ra"))
        & np.isfinite(as_float(joined, "secondary_dec"))
    )
    if np.any(coord_mask):
        q = SkyCoord(as_float(joined, "q_ra")[coord_mask] * u.deg, as_float(joined, "q_dec")[coord_mask] * u.deg)
        primary = SkyCoord(as_float(joined, "primary_ra")[coord_mask] * u.deg, as_float(joined, "primary_dec")[coord_mask] * u.deg)
        secondary = SkyCoord(as_float(joined, "secondary_ra")[coord_mask] * u.deg, as_float(joined, "secondary_dec")[coord_mask] * u.deg)
        old_primary_is_nearest[coord_mask] = primary.separation(q).arcsec <= secondary.separation(q).arcsec
        nearest_ok[coord_mask] = True

    old_primary_is_qso = np.where(can_use_colour, old_primary_is_qso_colour, old_primary_is_nearest)
    fallback_to_rband = ~can_use_colour & ~nearest_ok
    old_primary_is_qso[fallback_to_rband] = True

    label_source = np.full(len(joined), "qso_color", dtype=object)
    label_source[~can_use_colour & nearest_ok] = "qso_nearest_fallback"
    label_source[fallback_to_rband] = "rband_primary_fallback"

    out = joined.copy()
    out["qso_color_primary_is_rband_primary"] = old_primary_is_qso
    out["qso_color_label_source"] = label_source
    out["qso_color_reference_g_minus_r"] = med
    out["qso_color_reference_sigma_g_minus_r"] = sig
    out["qso_color_reference_p16_g_minus_r"] = p16
    out["qso_color_reference_p84_g_minus_r"] = p84
    out["qso_color_reference_n"] = n_ref
    out["qso_color_primary_distance_sigma"] = np.where(old_primary_is_qso, primary_scaled, secondary_scaled)
    out["qso_color_companion_distance_sigma"] = np.where(old_primary_is_qso, secondary_scaled, primary_scaled)
    out["rband_delta_g_minus_r"] = np.asarray(joined["delta_g_minus_r"], dtype=float)
    out["rband_delta_r_mag"] = np.asarray(joined["delta_r_mag"], dtype=float)

    swap_mask = ~old_primary_is_qso
    for first, second in (
        ("primary_ra", "secondary_ra"),
        ("primary_dec", "secondary_dec"),
        ("primary_gaia_source_id", "secondary_gaia_source_id"),
        ("primary_gaia_sep_arcsec", "secondary_gaia_sep_arcsec"),
        ("primary_gaia_g", "secondary_gaia_g"),
        ("primary_pm_sig", "secondary_pm_sig"),
        ("primary_parallax_over_error", "secondary_parallax_over_error"),
        ("primary_ruwe", "secondary_ruwe"),
        ("primary_is_original_gaia", "secondary_is_original_gaia"),
    ):
        swap_columns(out, first, second, swap_mask)

    old_delta_gr = np.asarray(joined["delta_g_minus_r"], dtype=float)
    old_delta_r = np.asarray(joined["delta_r_mag"], dtype=float)
    old_flux_ratio = np.asarray(joined["secondary_to_primary_flux_r"], dtype=float)
    out["delta_g_minus_r"] = np.where(old_primary_is_qso, old_delta_gr, -old_delta_gr)
    out["delta_r_mag"] = np.where(old_primary_is_qso, old_delta_r, -old_delta_r)
    new_flux_ratio = np.full(len(out), np.nan, dtype=float)
    valid_ratio = np.isfinite(old_flux_ratio) & (old_flux_ratio > 0.0)
    new_flux_ratio[valid_ratio & old_primary_is_qso] = old_flux_ratio[valid_ratio & old_primary_is_qso]
    new_flux_ratio[valid_ratio & ~old_primary_is_qso] = 1.0 / old_flux_ratio[valid_ratio & ~old_primary_is_qso]
    out["secondary_to_primary_flux_r"] = new_flux_ratio

    primary_sid = finite_source_id(out, "primary_gaia_source_id")
    secondary_sid = finite_source_id(out, "secondary_gaia_source_id")
    primary_has_gaia = primary_sid > 0
    secondary_has_gaia = secondary_sid > 0
    out["primary_has_gaia"] = primary_has_gaia
    out["secondary_has_gaia"] = secondary_has_gaia
    out["secondary_has_distinct_gaia"] = primary_has_gaia & secondary_has_gaia & (secondary_sid != primary_sid)

    summary_dict = {
        "method": "QSO-colour primary relabelling using redshift-local median g-r from nearest Tractor components.",
        "z_tol": float(z_tol),
        "n_input_rows": int(len(summary)),
        "n_qso_colour_labelled": int(np.sum(can_use_colour)),
        "n_qso_nearest_fallback": int(np.sum(~can_use_colour & nearest_ok)),
        "n_rband_primary_fallback": int(np.sum(fallback_to_rband)),
        "n_old_primary_retained": int(np.sum(old_primary_is_qso)),
        "n_old_secondary_promoted": int(np.sum(~old_primary_is_qso)),
        "median_reference_n": float(np.nanmedian(n_ref[np.isfinite(med)])) if np.any(np.isfinite(med)) else math.nan,
    }
    return out, summary_dict
