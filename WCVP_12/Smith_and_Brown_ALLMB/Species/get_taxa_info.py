import os

from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_columns, get_all_taxa

from library_info import FAMILIES_IN_GENTIANALES, WCVP_VERSION

_inputs_path = resource_filename(__name__, 'inputs')


def main():
    all_taxa = get_all_taxa(families_of_interest=FAMILIES_IN_GENTIANALES, accepted=True, version=WCVP_VERSION)

    species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']

    species_list = species_df[wcvp_columns['name']].unique().tolist()
    with open(os.path.join('inputs', 'species_list.txt'), 'w') as f:
        for line in species_list:
            f.write(f"{line}\n")

    species_df[[wcvp_columns['name'], wcvp_columns['family']]].to_csv(os.path.join('inputs', 'species_family_list.csv'))

if __name__ == '__main__':
    main()
