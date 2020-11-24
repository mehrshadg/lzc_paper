from neuro_helper.abstract.map import Space
from neuro_helper.map import ColeTemplateMap, T1T2Topo

__all__ = ["FMRI_RAW_DATA_ROOT_DIR", "cole_tpt", "cole_tpt_32", "t1t2_topo"]

FMRI_RAW_DATA_ROOT_DIR = "/data/hcp"
cole_tpt = ColeTemplateMap(Space.K59)()
cole_tpt_cortex = ColeTemplateMap(Space.K59_CORTEX)()
cole_tpt_32 = ColeTemplateMap(Space.K32)()
t1t2_topo = T1T2Topo(cole_tpt_cortex)()
