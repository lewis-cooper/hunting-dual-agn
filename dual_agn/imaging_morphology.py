"""Measure the original pair-axis/perpendicular r-band profile width ratio.

R_perp = sigma_along / sigma_perpendicular, using positive baseline-subtracted
profile moments. No visual decisions or per-object profile plots are involved.
"""

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


def axis_profile(image, along, perpendicular, lower, upper, bin_width, half_width=1.):
    edges = np.arange(lower, upper + .5 * bin_width, bin_width)
    if len(edges) < 2:
        raise ValueError('Profile range needs at least one bin')
    centers = .5 * (edges[:-1] + edges[1:])
    means = np.full(len(centers), np.nan)
    strip = np.isfinite(image) & (np.abs(perpendicular) <= half_width)
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        selected = strip & (along >= left) & (along < right)
        if selected.any():
            means[index] = image[selected].mean()
    return centers, means


def profile_sigma(coordinates, values):
    finite = np.isfinite(coordinates) & np.isfinite(values)
    coordinates, values = coordinates[finite], values[finite]
    if len(coordinates) < 3:
        return np.nan
    edge_width = .2 * max(np.ptp(coordinates), 1e-6)
    edges = (coordinates <= coordinates.min() + edge_width) | (coordinates >= coordinates.max() - edge_width)
    baseline_values = values[edges] if edges.sum() >= 3 else values
    baseline = np.median(baseline_values)
    weights = np.maximum(values - baseline, 0.)
    if (weights > 0).sum() < 2 or weights.sum() <= 0:
        return np.nan
    centroid = np.sum(weights * coordinates) / weights.sum()
    return float(np.sqrt(np.sum(weights * (coordinates - centroid)**2) / weights.sum()))


def measure_image(image, header, primary_ra, primary_dec, secondary_ra, secondary_dec):
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError('Morphology images require a celestial WCS')
    coordinates = SkyCoord([primary_ra, secondary_ra], [primary_dec, secondary_dec], unit='deg')
    component_x, component_y = wcs.world_to_pixel(coordinates)
    scale = float(header.get('PIXSCALE', 3600 * np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Image pixel scale must be finite and positive')
    shape_y, shape_x = image.shape
    inside = ((component_x >= 0) & (component_x < shape_x) &
              (component_y >= 0) & (component_y < shape_y))
    if not inside.all():
        raise ValueError('Both components must lie inside the supplied image')
    displacement = np.array([np.diff(component_x)[0], np.diff(component_y)[0]])
    distance_pixels = np.linalg.norm(displacement)
    if not np.isfinite(distance_pixels) or distance_pixels <= 0:
        raise ValueError('Morphology requires distinct finite component positions')
    direction_x, direction_y = displacement / distance_pixels
    pixel_y, pixel_x = np.indices(image.shape, dtype=float)
    relative_x = pixel_x - component_x[0]
    relative_y = pixel_y - component_y[0]
    along = (relative_x * direction_x + relative_y * direction_y) * scale
    perpendicular = (-relative_x * direction_y + relative_y * direction_x) * scale
    separation_arcsec = distance_pixels * scale
    rectangles = [(-1.5, separation_arcsec + 1.5, -1., 1.),
                  (.5 * separation_arcsec - 1., .5 * separation_arcsec + 1., -1.5, 1.5)]
    for lower_along, upper_along, lower_perp, upper_perp in rectangles:
        for axis_position in (lower_along, upper_along):
            for perpendicular_position in (lower_perp, upper_perp):
                corner_x = component_x[0] + (axis_position * direction_x - perpendicular_position * direction_y) / scale
                corner_y = component_y[0] + (axis_position * direction_y + perpendicular_position * direction_x) / scale
                if not (-.5 <= corner_x <= shape_x - .5 and -.5 <= corner_y <= shape_y - .5):
                    raise ValueError('Image does not cover the full morphology profile apertures')
    aperture = (((along >= -1.5) & (along <= separation_arcsec + 1.5) & (np.abs(perpendicular) <= 1.)) |
                ((np.abs(along - .5 * separation_arcsec) <= 1.) & (np.abs(perpendicular) <= 1.5)))
    if not np.isfinite(image[aperture]).all():
        raise ValueError('Image has missing pixels inside the morphology profile apertures')
    axis = axis_profile(image, along, perpendicular, -1.5, separation_arcsec + 1.5, scale)
    control = axis_profile(image, perpendicular, along - .5 * separation_arcsec, -1.5, 1.5, scale)
    sigma_axis, sigma_perpendicular = profile_sigma(*axis), profile_sigma(*control)
    ratio = sigma_axis / sigma_perpendicular if np.isfinite(sigma_perpendicular) and sigma_perpendicular > 0 else np.nan
    return dict(image_axis_sigma_arcsec=sigma_axis, image_perp_sigma_arcsec=sigma_perpendicular,
                r_perp_sigma_ratio=ratio, morphology_pixel_scale_arcsec=scale,
                morphology_separation_arcsec=separation_arcsec)
