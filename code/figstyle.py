"""Shared figure style for the LDG project."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'  # keep text editable in SVG
matplotlib.rcParams['figure.dpi'] = 120
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['font.size'] = 9

# colorblind-friendly habitat palette (Okabe-Ito inspired, extended)
HABITAT_COLORS = {
    'soil_natural':       '#009E73',
    'soil_agricultural':  '#56B4E9',
    'marine_planktonic':  '#0072B2',
    'marine_sediment':    '#332288',
    'freshwater':         '#88CCEE',
    'wastewater':         '#CC6677',
    'urban_built':        '#DDCC77',
    'plant_rhizosphere':  '#44AA99',
    'plant_phyllosphere': '#999933',
    'animal_gut':         '#E69F00',
    'animal_skin':        '#F0E442',
    'human_gut':          '#D55E00',
}

HABITAT_LABELS = {
    'soil_natural': 'Soil (natural)',
    'soil_agricultural': 'Soil (agricultural)',
    'marine_planktonic': 'Marine (planktonic)',
    'marine_sediment': 'Marine (sediment)',
    'freshwater': 'Freshwater',
    'wastewater': 'Wastewater',
    'urban_built': 'Urban / built',
    'plant_rhizosphere': 'Plant rhizosphere',
    'plant_phyllosphere': 'Plant phyllosphere',
    'animal_gut': 'Animal gut',
    'animal_skin': 'Animal skin/surface',
    'human_gut': 'Human gut',
}

# climate exposure scores (plan section 4.4): free-living=3, plant=2, host=1, engineered=0
EXPOSURE = {
    'soil_natural': 3, 'soil_agricultural': 3, 'marine_planktonic': 3,
    'marine_sediment': 3, 'freshwater': 3,
    'plant_rhizosphere': 2, 'plant_phyllosphere': 2,
    'animal_gut': 1, 'animal_skin': 1, 'human_gut': 1,
    'wastewater': 0, 'urban_built': 0,
}
