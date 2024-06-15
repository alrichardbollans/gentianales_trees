import os

from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_columns, get_all_taxa

from library_info import FAMILIES_IN_GENTIANALES, WCVP_VERSION, MY_MAIN_FAMILIES

_inputs_path = resource_filename(__name__, 'inputs')


def main():
    all_taxa = get_all_taxa(families_of_interest=FAMILIES_IN_GENTIANALES, accepted=True, version=WCVP_VERSION)

    genus_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Genus']

    genus_list = genus_df[wcvp_columns['name']].unique().tolist()
    with open(os.path.join('inputs', 'genus_list.txt'), 'w') as f:
        for line in genus_list:
            f.write(f"{line}\n")

    genus_df[[wcvp_columns['name'], wcvp_columns['family']]].to_csv(os.path.join('inputs', 'genus_family_list.csv'))


    main_family_df = genus_df[genus_df[wcvp_columns['family']].isin(MY_MAIN_FAMILIES)]
    genus_list = main_family_df[wcvp_columns['name']].unique().tolist()
    with open(os.path.join('inputs', 'apoc_log_rub_genus_list.txt'), 'w') as f:
        for line in genus_list:
            f.write(f"{line}\n")

    main_family_df[[wcvp_columns['name'], wcvp_columns['family']]].to_csv(os.path.join('inputs', 'apoc_log_rub_genus_family_list.csv'))


if __name__ == '__main__':
    main()
