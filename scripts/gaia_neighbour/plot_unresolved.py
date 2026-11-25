import numpy as np
import pandas as pd
from dual_agn.plot_style import plt, save_figure
from dual_agn.settings import IDENTIFIER_DTYPES, ONE_ROOT

systems = pd.read_csv(
    ONE_ROOT / 'scoring_inputs.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)
blend_ranked = pd.read_csv(
    ONE_ROOT / 'ranked_unresolved.csv',
    dtype=IDENTIFIER_DTYPES,
    float_precision='round_trip'
)

fig = plt.figure(figsize=(11, 6))
ax = fig.add_subplot(111, projection='aitoff')
colours = plt.get_cmap('tab10').colors
for index, (survey, survey_systems) in enumerate(systems.groupby('primary_qso_catalog', sort=True)):
    ra_wrapped_deg = (survey_systems['primary_qso_ra'] + 180) % 360 - 180
    ax.scatter(
        -np.deg2rad(ra_wrapped_deg),
        np.deg2rad(survey_systems['primary_qso_dec']),
        s=18,
        alpha=0.7,
        color=colours[index % len(colours)],
        label=f'{survey} ({len(survey_systems)})'
    )
ra_tick_longitudes = np.arange(-120, 121, 60)
ax.set_xticks(np.deg2rad(ra_tick_longitudes))
ax.set_xticklabels([f'{-longitude % 360}°' for longitude in ra_tick_longitudes])
ax.set_yticks(np.deg2rad(np.arange(-60, 61, 30)))
ax.set(xlabel='RA', ylabel='Dec')
ax.set_title(f'Sky positions of {len(systems)} candidate systems', pad=25)
ax.grid(alpha=0.3)
ax.legend(title='Source', loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3, frameon=False)
save_figure(fig, 'input_sky', 'one_arcsecond')
plt.close("all")


unscored_systems = blend_ranked.loc[
    ~np.isfinite(blend_ranked['log_bf_qq_qs_magdep']),
    [
        'qso_group_id',
        'primary_qso_catalog',
        'primary_qso_z',
        'status_magdep',
        'decals_system_assignment',
        'blend_r_snr',
        'fallback_r_band_input',
        'fallback_r_snr',
        'galactic_abs_b_deg'
    ]
].merge(
    systems[['qso_group_id', 'primary_qso_ra', 'primary_qso_dec']],
    on='qso_group_id',
    validate='one_to_one'
)
unscored_systems['redshift_source'] = unscored_systems['primary_qso_catalog'].replace({'DESI_DR1_QSO': 'DESI', 'SDSS_DR16Q': 'SDSS', 'QUAIA_G20_V1': 'Quaia'})
legacy_r = unscored_systems['decals_system_assignment'].map({
    'primary': 'Legacy r (primary match)',
    'secondary': 'Legacy r (secondary match)',
    'desi_target': 'Legacy r (DESI target)'
})
fallback_r = unscored_systems['fallback_r_band_input'].map({'ps1:r': 'PS1 r', 'nsc:r': 'NSC r', 'skymapper:r': 'SkyMapper r'})
unscored_systems['r_band_for_scoring'] = legacy_r.fillna(fallback_r).fillna('none')
unscored_systems['r_snr'] = unscored_systems['blend_r_snr'].where(legacy_r.notna(), unscored_systems['fallback_r_snr'])
unscored_systems['why_unscorable'] = unscored_systems['status_magdep'].map({
    'outside_latitude': '|b| < 25°',
    'no_reference_band': 'No usable r measurement',
    'redshift_outside_support': 'Redshift outside model support',
    'empty_hypotheses': 'No viable QQ or QS model hypothesis',
    'too_faint_for_two_components': 'r S/N < 20'
}).fillna(unscored_systems['status_magdep'])
low_latitude = unscored_systems['galactic_abs_b_deg'].lt(25)
no_r = unscored_systems['r_band_for_scoring'].eq('none')
weak_r = unscored_systems['r_snr'].lt(20)
unscored_systems.loc[low_latitude & no_r, 'why_unscorable'] += '; no usable r'
unscored_systems.loc[low_latitude & weak_r, 'why_unscorable'] += '; r S/N < 20'
unscored_systems = unscored_systems.rename(columns={
    'primary_qso_ra': 'ra_deg',
    'primary_qso_dec': 'dec_deg',
    'primary_qso_z': 'redshift'
})[[
    'qso_group_id',
    'ra_deg',
    'dec_deg',
    'redshift',
    'redshift_source',
    'r_band_for_scoring',
    'r_snr',
    'galactic_abs_b_deg',
    'why_unscorable'
]].sort_values(['why_unscorable', 'qso_group_id'])
print(f'Unscored systems: {len(unscored_systems)}')
unscored_systems.to_csv(ONE_ROOT / 'unscored_unresolved_diagnostic.csv', index=False)
