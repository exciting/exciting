""" Module to change GROUP settings in  defaults_config.yml,
typically used when running an instance in the CI.

Safer than using sed because one can easily throw an error
 - "sed -i 's/  NONE.*/  NONE: False/g' defaults_config.yml" 
 - "sed -i 's/  GW.*/  GW: True/g' defaults_config.yml"
"""
import os
import yaml 
import argparse as ap
from typing import Tuple


def parse_key_changes() -> Tuple[str, dict]:
    """ Parse keys specified at command line.
    """
    p = ap.ArgumentParser(
        description="Expects -file and -groups. "
        "For example: python3 set_groups.py -file defaults_config.yml -groups GW BSE")
    
    p.add_argument("-groups",type=str, nargs='+', help='Groups to run', default = [], required=False)

    p.add_argument('-file',type=str, help='Defaults YAML file', required=True)

    groups_to_run = {}
    for group in p.parse_args().groups:
        groups_to_run[group] = True

    return p.parse_args().file, groups_to_run


def modify_group_settings(file_name, groups_to_run: dict):
    """ Modify default group settings in the YAML config file.
    """
    if not os.path.exists(file_name):
        raise FileExistsError(f'{file_name} cannot be found')

    with open(file_name, "r") as stream:
        try:
            defaults = yaml.safe_load(stream)
        except yaml.YAMLError:
            raise yaml.YAMLError()
        
    # Set all groups to not run
    group_settings = {name: False for name in defaults['group_execution'].keys()}

    # Modify groups to run
    for group in groups_to_run:
        if group not in group_settings:
            raise KeyError(f'{group} is not a valid group in {file_name}. '
                           f'Valid groups are: {list(group_settings)}')
        group_settings[group] = True

    defaults['group_execution'] = group_settings
    
    with open(file_name, 'w') as fid:
        yaml.dump(defaults, fid, default_flow_style=False)
    
    print(f'Running groups: {groups_to_run}')


if __name__ == "__main__":
    file_name, groups_to_run = parse_key_changes()
    modify_group_settings(file_name, groups_to_run)
