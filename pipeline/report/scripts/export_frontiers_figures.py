"""Export manuscript figures for Frontiers: 180 mm wide, 300 dpi TIFF, in order of appearance."""
from pathlib import Path
from PIL import Image

REPORT = Path(__file__).resolve().parents[1]
FIGURES = REPORT.parent / 'outputs' / 'figures'
OUT = REPORT / 'submission' / 'figures'
ORDER = ['fig35_manuscript_gbr_and_workflow', 'fig23_workflow_decision', 'fig02_four_reef_comparison',
         'fig13_manta_aggregation_cv', 'fig20_high_density_retention', 'fig34_ghost_per_miss',
         'fig22_metacommunity_cots_coral', 'fig24_metacommunity_additional_reefs',
         'fig27_metacommunity_projection_effects']
DPI = 300
MAX_WIDTH = round(180 / 25.4 * DPI)

OUT.mkdir(parents=True, exist_ok=True)
for number, name in enumerate(ORDER, 1):
    image = Image.open(FIGURES / f'{name}.png').convert('RGB')
    if image.width > MAX_WIDTH:
        image = image.resize((MAX_WIDTH, round(image.height * MAX_WIDTH / image.width)), Image.LANCZOS)
    path = OUT / f'Figure_{number}.tif'
    image.save(path, dpi=(DPI, DPI), compression='tiff_lzw')
    print(path.name, image.size)
