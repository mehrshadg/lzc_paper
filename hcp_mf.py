import numpy as np
import glob

from neuro_helper.abstract.map import TemplateMap
from neuro_helper.measurement import calc_mf
from joblib import Parallel, delayed
from neuro_helper.hcp.fmri.storage import load_raw_files, FMRILocalStorage, concat_scans
from neuro_helper.hcp.fmri.generic import task_order
from neuro_helper.statistics import fir_filter, welch_psd
from neuro_helper.generic import out_of, generate_long_data
from neuro_helper.generic import find_shared_subjects as fs_subjects
from neuro_helper.storage import ANYTHING
import config


def do_a_subject(tpt: TemplateMap, input):
    data, fs, label = load_raw_files(input, tpt.space, concat_scans)
    print(f"Calculating MF on {label}")

    data, freq_l, freq_h = fir_filter(data, fs, max_freq_low=0.01, pass_type="hp")
    data = tpt.parcellate(data)
    freq, psd = welch_psd(data, freq_l, fs)
    filt = freq >= freq_l
    if freq_h is not None:
        filt = filt & freq <= freq_h
    return calc_mf(freq[filt], psd[:, filt])


def run_script(tpt: TemplateMap):
    for task in task_order(tpt.space, with_rest=True):
        storage = FMRILocalStorage(config.FMRI_RAW_DATA_ROOT_DIR, task, ANYTHING, tpt.space)
        files_dict = storage.get_all_by_scan()
        for key, file_infos in files_dict.items():
            output_file = out_of(f"fmri-hcp-{task}.mf.rois-{tpt.name}.scan-{key}-{tpt.space}.npy", False)
            subj_ids = list(file_infos.keys())
            output = np.asarray(Parallel(n_jobs=30)(delayed(do_a_subject)(
                tpt, file_infos[subj_id]) for subj_id in subj_ids))
            np.save(output_file, (task, key, subj_ids, output))


def find_files(**kwargs):
    task = kwargs["task"]
    tpt = kwargs["template"]
    files = glob.glob(out_of(f"fmri-hcp-{task}.mf.rois-{tpt.name}.scan-{ANYTHING}-{tpt.space}.npy", False))
    files.sort()
    return files


def prepare_file_content(ret_metric):
    if ret_metric:
        return lambda content: (content[1], content[2], content[3][:, 1:])
    else:
        return lambda content: (content[1], content[2])


def find_shared_subjects(tpt: TemplateMap, tasks, return_indices=False):
    return fs_subjects(find_files, prepare_file_content(False), tpt, tasks, return_indices)


def gen_long_data(tpt: TemplateMap, show_warning=False):
    return generate_long_data(find_files, prepare_file_content(True), tpt, task_order(tpt.space, True), show_warning)


if __name__ == "__main__":
    run_script(config.cole_tpt)
