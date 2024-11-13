import os

import numpy as np
import pandas as pd
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_columns, get_all_taxa, wcvp_accepted_columns
from wcvpy.wcvp_name_matching import get_genus_from_full_name

from WCVP_12.v12_library_info import FAMILIES_IN_GENTIANALES, WCVP_VERSION

_inputs_path = resource_filename(__name__, 'inputs')


# TODO: Use name matching to get species.relative information
def make_hybrid_epithet_one_word(name):
    ## Doesn't seem like uphylomaker can handle hybrids
    if name is not None and name==name:
        if name.startswith('×'):
            return np.nan
            # return name.replace('× ', '×_')
        else:
            return name
    else:
        return name

def main():
    # Get sp relative information using synonym information from original smb tree
    name_matching_file_from_smb = pd.read_csv(os.path.join('..', 'Smith_and_Brown_ALLMB','Genus', 'inputs', 'acc_name_tree_Gentianales.csv'))

    name_matching_file_from_smb = name_matching_file_from_smb[name_matching_file_from_smb['taxon_status'] == 'Synonym']
    name_matching_file_from_smb = name_matching_file_from_smb[name_matching_file_from_smb['accepted_rank'] == 'Species']
    name_matching_file_from_smb = name_matching_file_from_smb[['binomial_name', 'accepted_species']]

    # Some accepted names that aren't in the tree are resolved to by multiple tree synonyms
    # but only one can be used by U phylomaker
    name_matching_file_from_smb = name_matching_file_from_smb.drop_duplicates(subset=['accepted_species'])

    name_matching_file_from_smb = name_matching_file_from_smb.rename(columns={'binomial_name': 'species.relative',
                                                                              'accepted_species': 'species'})
    name_matching_file_from_smb['genus.relative'] = name_matching_file_from_smb['species.relative'].apply(get_genus_from_full_name)
    all_taxa = get_all_taxa(families_of_interest=FAMILIES_IN_GENTIANALES, accepted=True, version=WCVP_VERSION)

    species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']

    species_df = species_df[[wcvp_accepted_columns['name'], wcvp_accepted_columns['parent_name'],
                             wcvp_accepted_columns['family']]]

    species_df = species_df.rename(columns={wcvp_accepted_columns['name']: 'species', wcvp_accepted_columns['parent_name']: 'genus',
                                            wcvp_accepted_columns['family']: 'family'})
    species_df = species_df.reset_index(drop=True)

    species_df = pd.merge(species_df, name_matching_file_from_smb, on='species', how='left')
    species_df['species'] = species_df['species'].apply(make_hybrid_epithet_one_word)
    species_df = species_df.dropna(subset=['species'])
    species_df.to_csv(os.path.join('inputs', 'species_family_list.csv'), index=False)

    genus_df = species_df[['genus', 'family']]
    genus_df = genus_df.drop_duplicates(keep='first')
    genus_df = genus_df.reset_index(drop=True)

    genus_df.to_csv(os.path.join('inputs', 'genus_family_list.csv'), index=False)


if __name__ == '__main__':
    main()
