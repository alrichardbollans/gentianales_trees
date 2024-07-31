import os

from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_columns, get_all_taxa, wcvp_accepted_columns

from WCVP_12.v12_library_info import FAMILIES_IN_GENTIANALES, WCVP_VERSION

_inputs_path = resource_filename(__name__, 'inputs')


def main():
    all_taxa = get_all_taxa(families_of_interest=FAMILIES_IN_GENTIANALES, accepted=True, version=WCVP_VERSION)

    species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']

    species_df = species_df[[wcvp_accepted_columns['name'], wcvp_accepted_columns['parent_name'],
                             wcvp_accepted_columns['family']]]

    species_df = species_df.rename(columns={wcvp_accepted_columns['name']: 'species', wcvp_accepted_columns['parent_name']: 'genus',
                                            wcvp_accepted_columns['family']: 'family'})
    species_df = species_df.reset_index(drop=True)
    species_df.to_csv(os.path.join('inputs', 'species_family_list.csv'), index=False)

    genus_df = species_df[['genus', 'family']]
    genus_df = genus_df.drop_duplicates(keep='first')
    genus_df = genus_df.reset_index(drop=True)

    genus_df.to_csv(os.path.join('inputs', 'genus_family_list.csv'), index=False)


if __name__ == '__main__':
    main()
