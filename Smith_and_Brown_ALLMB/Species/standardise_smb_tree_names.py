import os
import re

import pandas as pd
from library_info import WCVP_VERSION
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_accepted_columns
from wcvpy.wcvp_name_matching import get_accepted_info_from_names_in_column

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')

_output_path = resource_filename(__name__, 'outputs')

tree_file = os.path.join('..', 'SMB_ALLMB_Gentianales.tre')  # SMB tree from gentianales node.


def get_binomial_from_label(leaf: str):
    """

    This method extracts the binomial name from a given label string. The label may contain additional information such as epithet indicators or other substrings. The binomial name consists of two parts: the genus and the specific epithet.

    :param leaf: The label from which the binomial name will be extracted.
    :return: The extracted binomial name.

    Example usage:
    ```python
    leaf_label = "Canis_lupus_subsp."
    binomial_name = get_binomial_from_label(leaf_label)
    print(binomial_name)
    # Output: "Canis lupus"
    ```

    Note that any substring after the epithet indicators will be removed, and multiple words in the label will be joined using spaces.

    """
    # Remove everything after these
    epithet_indicators = ['_subsp.', '_sp.', '_var.']
    for e in epithet_indicators:
        leaf = leaf.split(e, 1)[0]
    x = " ".join(leaf.split(sep='_')[:2])
    return x.strip()


def substitute_name_in_tree(tree_string: str, old_name: str, new_name: str):
    pattern = r'(?<=[(),])({})(?=;|:)'.format(re.escape(old_name))

    # Note ape will not read spaces, so add underscores back to names
    tree_string = re.sub(pattern, new_name.replace(' ', '_'), tree_string)

    return tree_string


def relabel_tree(families_of_interest: list, outfile: str):
    f = open(tree_file, "r")
    tree_string = f.readline()
    # From https://stackoverflow.com/questions/45668107/python-regex-parsing-newick-format
    rx = r'[(),]+([^;:]+)(?=;|:)'  # Looks for any combination of ((,, greedily, then ahead matches any non ;: characters and terminates with a ; or :
    name_list = re.findall(rx, tree_string)
    binomial_names = [get_binomial_from_label(x) for x in name_list]
    zipped = list(zip(name_list, binomial_names))

    df = pd.DataFrame(zipped, columns=['tree_name', 'binomial_name'])

    acc_name_df = get_accepted_info_from_names_in_column(df, 'binomial_name', wcvp_version=WCVP_VERSION)
    acc_name_df.to_csv(os.path.join('inputs', 'acc_name_tree_Gentianales.csv'))
    acc_name_df = pd.read_csv(os.path.join('inputs', 'acc_name_tree_Gentianales.csv'), index_col=0)

    # Catch words in tree string by left hand word boundaries ((),) and right hand ; or : characters
    for index, row in acc_name_df.iterrows():
        print(f'{index} out of {len(acc_name_df)}')
        if row[wcvp_accepted_columns['family']] in families_of_interest and row['accepted_rank'] == 'Species':
            # Note ape will not read spaces, so add underscores back to names
            tree_string = substitute_name_in_tree(tree_string, row['tree_name'], row[wcvp_accepted_columns['name']])


        elif row['tree_name'] != 'Gentianales.rn.d8s.tre':
            tree_string = substitute_name_in_tree(tree_string, row['tree_name'], 'NON_FAMILY_TIP')

    f = open(outfile, "w")
    f.writelines([tree_string])


def main():
    standard_apoc_log_rub_file = os.path.join(_temp_outputs_path, 'standardised_apoc_log_rub_smb_tree.tre')

    relabel_tree(['Apocynaceae', 'Loganiaceae', 'Rubiaceae'], standard_apoc_log_rub_file)


if __name__ == '__main__':
    main()
