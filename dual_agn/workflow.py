"""Pipeline choices and their public product files."""

PIPELINE_CHOICES = ('both', 'gaia_neighbour', 'imaging')


def selected_pipelines(pipeline):
    if pipeline not in PIPELINE_CHOICES:
        raise ValueError('Choose pipeline both, gaia_neighbour or imaging')
    return {'gaia_neighbour', 'imaging'} if pipeline == 'both' else {pipeline}


def product_names(pipeline):
    names = {'final_candidates.csv', 'threshold_calibration.json'}
    if 'gaia_neighbour' in selected_pipelines(pipeline):
        names.update({'one_arcsecond_scores.csv', 'one_arcsecond_cutflow.csv'})
    if 'imaging' in selected_pipelines(pipeline):
        names.update({'imaging_scores.csv', 'imaging_cutflow.csv'})
    return names
