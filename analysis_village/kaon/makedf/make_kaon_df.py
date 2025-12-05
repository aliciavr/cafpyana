"""
Kaon analysis data frame maker
"""

import numpy as np
import pandas as pd

import makedf.makedf as makedf
import makedf.branches as branches
import pyanalib.pandas_helpers as ph
import makedf.util as util


KPDG = {
    'kplus': makedf.PDG["kaon_p"][0],
    'kzero': makedf.PDG["kaon_0"][0]
}
KMASS = {
    'kplus': makedf.PDG["kaon_p"][2],
    'kzero': makedf.PDG["kaon_0"][2]
}

# TODO pick a better number
TRUE_KE_CUT = 0.

# add some extra branches into primary ones, use set to make sure it stays
# unique even if mcprimbranches changes in the future
PRIM_BRANCHES = list(set(branches.mcprimbranches + [
    'rec.mc.nu.prim.G4ID',
    'rec.mc.nu.prim.end_process',
]))


def signal(df: pd.DataFrame, ktype: str, cc: bool=True) -> pd.DataFrame:
    # InFV requires inzback but it does nothing for SBND case
    # put NaN here so we'll hopefully get an error if this ever changes
    is_true_fv = util.InFV(df.position, inzback=np.nan, det='SBND')

    # TODO daughter particle containment cut

    has_k = (getattr(df, f'n{ktype}') > 0)
    cc_nc = (df.iscc == cc)
    return is_true_fv & cc_nc & has_k


def make_kaon_mcdf(f: pd.DataFrame) -> pd.DataFrame:
    mcdf = makedf.make_mcdf(f)

    mcprimdf = ph.loadbranches(f["recTree"], PRIM_BRANCHES).rec.mc.nu.prim

    # for daughter info, will match on entry + G4ID
    tpartdf = ph.loadbranches(f["recTree"], branches.trueparticlebranches).rec.true_particles
    tpartdf = tpartdf.reset_index().set_index(['entry', 'G4ID'])

    # add kaon info: Number above threshold & primary info
    for kname in ('kplus', 'kzero'):
        # number of kaons above KE threshold
        ke = mcprimdf[mcprimdf.pdg==KPDG[kname]].genE - KMASS[kname] 
        mcdf = ph.multicol_add(mcdf, ((mcprimdf.pdg==KPDG[kname]) \
                                              & (ke > TRUE_KE_CUT)).groupby(level=[0, 1]).sum().rename(f'n{kname}'))

        # kaon primary info
        kdf = mcprimdf[mcprimdf.pdg==KPDG[kname]]


        # kaon daughter info
        mcprimdaughtersdf = makedf.make_mcprimdaughtersdf(f).rec.mc.nu.prim
        kdaughterdf = mcprimdaughtersdf[mcprimdaughtersdf.index.droplevel(-1).isin(kdf.index)]
        kdaughter_tpartdf = tpartdf[
            tpartdf.index.isin(pd.MultiIndex.from_frame(kdaughterdf.reset_index()[['entry', 'daughters']]))
        ].reset_index().set_index(['entry', 'parent'])

        print(kdaughter_tpartdf.groupby(level=[0, 1]).size().rename('ndaughters'))

        kdf.columns = pd.MultiIndex.from_tuples([tuple([kname] + list(c)) for c in kdf.columns])
        mcdf = ph.multicol_merge(mcdf, kdf, how="left", left_index=True, right_index=True, validate="one_to_one")

    mcdf['is_signal_kp_cc'] = signal(mcdf, ktype='kplus', cc=True)
    mcdf['is_signal_kp_nc'] = signal(mcdf, ktype='kplus', cc=False)

    return mcdf


def make_kaon_recodf(f: pd.DataFrame) -> pd.DataFrame:
    pandora_df = makedf.make_pandora_df(f)

    # precuts
    pandora_df = pandora_df[util.InFV(pandora_df.slc.vertex, inzback=np.nan, det='SBND')]
    pandora_df = pandora_df[pandora_df.slc.is_clear_cosmic == 0]

    return pandora_df
