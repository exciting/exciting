import os
import requests

from argparse import ArgumentParser
from io import BytesIO
from zipfile import ZipFile
from pathlib import Path

CALCULATION_ENTRY_MAP = {
    'groundstate': 'rdiYhF-9QnBlL3IDu3y0OV8ZXnuJ',
    'metaGGA': '_-3Rc4A3VfNeOvt3SadrxnB711g3',
    'GW': 'HUMNIOp4pf8ueAGjmviJp66Pqd4i',
    'phonons': '-CPyhiFDrYXaHta7io18x-_QfmLT',
    'BSE': 'Xr7py4FerQtXIXeK61j81UbQ2od8',
    'TDDFT': 'lyYhRaJM2KDxsxftrxBfYxlYc-hd',
    'hybrid': '_YwxC-ut5frhpQEPz7SaEq5MEBr4'
}

def get_entry_from_NOMAD(entry_id, base_url='https://nomad-lab.eu/prod/v1/api/v1', overwrite=False):
    """
    Get entry from NOMAD and extract to current directory.

    If the mainfile for this entry is not in the lowest directory, all intermediate folders will also be created.
    All additional folders created by NOMAD, as well as the metadata, are removed.
    """
    r = requests.get(f'{base_url}/entries/{entry_id}/raw')
    assert r.ok, f'Could not download data from NOMAD, failed with status code {r.status_code}: {r.reason}.'
    
    data_file_buffer = BytesIO(r.content)
    
    with ZipFile(data_file_buffer) as zip_ref:
        for zip_info in zip_ref.filelist:
            file_path = Path(zip_info.filename)
            if len(file_path.parent.name)>0:
                file_path = file_path.relative_to(file_path.parts[0])
                os.makedirs(file_path.parent, exist_ok=True)
                if os.path.exists(file_path) and not overwrite:
                    raise FileExistsError(f'{file_path} exists, set --overwrite flag to overwrite')
                with zip_ref.open(zip_info) as source, open(file_path, 'wb') as target:
                    target.write(source.read())

def main():
    parser = ArgumentParser(description='Download the data for for an exciting tutorial from NOMAD and extract it.')

    parser.add_argument('--base-url',
                        type=str,
                        default='https://nomad-lab.eu/prod/v1/api/v1',
                        help='URL of the NOMAD installation.')
    
    parser.add_argument('--overwrite',
                        action='store_true')
    
    parser.add_argument('CALC',
                        type=str,
                        help='Calculation type for which the data should be downloaded.',
                        choices=list(CALCULATION_ENTRY_MAP.keys()))
    
    args = parser.parse_args()


    get_entry_from_NOMAD(CALCULATION_ENTRY_MAP[args.CALC], args.base_url, overwrite=args.overwrite)
    
if __name__=='__main__':
    main()