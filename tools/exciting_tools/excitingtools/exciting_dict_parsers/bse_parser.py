"""Parsers for BSE output files.
"""
import re
from typing import Optional

import numpy as np


def numpy_gen_from_txt(name: str, skip_header: Optional[int] = 0) -> np.ndarray:
    """  Numpy genfromtxt, dressed in try/expect.

    Not worth generalising, as would need to support genfromtxt's API.

    :param name: File name.
    :param skip_header: Optional number of header lines to skip.
    :return data: Parsed data.
    """
    try:
        data = np.genfromtxt(name, skip_header=skip_header)
    except ValueError:
        raise ValueError(f'Failed to parse {name}')
    return data


def parse_EPSILON_NAR(name: str) -> dict:
    """
    Parser for:
        EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC.OUT.xml,
        EPSILON_NAR_FXCMB1_OC_QMT001.OUT.xml,
        EPSILON_NAR_NLF_FXCMB1_OC_QMT001.OUT.xml,
        LOSS_NAR_FXCMB1_OC_QMT001.OUT.xml
    """
    data = numpy_gen_from_txt(name, skip_header=14)
    out = {
        "frequency": data[:, 0],
        "real_oscillator_strength": data[:, 1],
        "imag_oscillator_strength": data[:, 2],
        "real_oscillator_strength_kkt": data[:, 3]
    }
    return out


def parse_LOSS_NAR(name):
    """
    Parser for:
     LOSS_NAR_FXCMB1_OC_QMT001.OUT.xml,
     LOSS_NAR_NLF_FXCMB1_OC_QMT001.OUT.xml
    """
    data = numpy_gen_from_txt(name, skip_header=14)
    out = {
        "frequency": data[:, 0],
        "real_oscillator_strength": data[:, 1],
        "imag_oscillator_strength": data[:, 2]
    }

    return out


def parse_EXCITON_NAR_BSE(name):
    """
    Parser for EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC.OUT
    """
    data = numpy_gen_from_txt(name, skip_header=14)
    out = {}
    out["state"] = data[:, 0]
    out["energy"] = data[:, 1]
    out["energy_shifted"] = data[:, 2]
    out["abs_oscillator_strength"] = data[:, 3]
    out["real_oscillator_strength"] = data[:, 4]
    out["imaginary_oscillator_strength"] = data[:, 5]

    return out


def parse_infoxs_out(name: str, parse_timing: bool = False) -> dict:
    """
    Parser for INFOXS.OUT file. Parses only the started and stopped tasks.
    Searches for lines like:
        'EXCITING <version> started for task <taskname> (<tasknumber>)'
    and
        'EXCITING <version> stopped for task <tasknumber>'
    See example file: exciting/test/test_farm/BSE/PBE_SOL-LiF/ref/INFOXS.OUT
    If a started task is found, it gets stored with name, number and status.
    If the task is found to be finished afterwards, the status finished is set to True.

    For success, the last started tasks has to be finished after that (in the file).
    Last finished task is the last task if calculation was successful, the task before that 
    if it finished, else None.
    :param name: path of the file to parse
    :param parse_timing: parse also timing information for the tasks. By default this is set to
                         False. If the task has not finished None is returned as timing.
    :returns: dictionary containing parsed file
    """
    with open(name) as file:
        lines = file.readlines()

    tasks = []
    current_task = -1

    lines = "\n".join(lines)
    all_tasks = re.findall(r'EXCITING .* (started) for task (.*) \( ?(\d+)\)|'
                           r'EXCITING .* stopped for task .* (\d+)', lines)

    for task in all_tasks:
        if task[0] == 'started':
            tasks.append({
                'name': task[1],
                'number': int(task[2]),
                'finished': False
            })
            current_task += 1
        else:
            # asserts shouldn't happen with Exciting:
            assert tasks != [], 'No tasks started!'
            assert tasks[current_task]['number'] == int(task[3]), 'Wrong task stopped.'
            tasks[current_task]['finished'] = True

    success = tasks[-1]['finished']
    last_finished_task = None
    if success:
        last_finished_task = tasks[-1]['name']
    elif len(tasks) > 1 and tasks[-2]['finished']:
        last_finished_task = tasks[-2]['name']

    if parse_timing:
        times = parse_times(lines)
        finished_tasks = [task for task in tasks if task['finished']]
        assert len(times['cpu']) == len(finished_tasks), 'Numbers of finished tasks and parsed times are not the same.'

        for index, task in enumerate(finished_tasks):
            task['cpu_time'] = float(times['cpu'][index])
            task['wall_time'] = float(times['wall'][index])
            task['cpu_time_cum'] = float(times['cpu_cum'][index])
            task['wall_time_cum'] = float(times['wall_cum'][index])
    
    return {'tasks': tasks,
            'success': success,
            'last_finished_task': last_finished_task}


def parse_times(infoxs_string: str) -> dict:
    """Parse the run times in INFOXS.OUT for each task.
    :param infoxs_string: String that contains the INFOXS.OUT file.
    :returns: dictionary containing a list of run times for each measurement.
    """
    cpu_times = re.findall(r'CPU time \s*: ([\d\.\d]+) sec', infoxs_string)
    wall_times = re.findall(r'wall time \s*: ([\d\.\d]+) sec', infoxs_string)
    cpu_times_cum = re.findall(r'CPU time \s* \(cumulative\) \s*: ([\d\.\d]+) sec', infoxs_string)
    wall_times_cum = re.findall(r'wall time \(cumulative\) \s*: ([\d\.\d]+) sec', infoxs_string)

    assert len(cpu_times) == len(wall_times), 'Numbers of parsed timings are not consistent.'
    assert len(cpu_times) == len(cpu_times_cum), 'Numbers of parsed timings are not consistent.'
    assert len(cpu_times) == len(wall_times_cum), 'Numbers of parsed timings are not consistent.'
    
    return {'cpu': cpu_times, 
            'wall': wall_times,
            'cpu_cum': cpu_times_cum,
            'wall_cum': wall_times_cum}
