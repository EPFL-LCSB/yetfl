# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:35:26 2020

To avoid using context specific coding in cobra function

@author: Omid
"""
import pandas as pd
from cobra.manipulation.delete import find_gene_knockout_reactions
from functools import partial
from itertools import product


def _reactions_knockouts_with_restore(model, reactions):
    prev_lb = dict()
    prev_ub = dict()
    for rxn in reactions:
        prev_lb[rxn.id] = rxn.lower_bound
        prev_ub[rxn.id] = rxn.upper_bound
        rxn.knock_out()
    growth = model.slim_optimize()
    for rxn in reactions:
        rxn.lower_bound = prev_lb[rxn.id]
        rxn.upper_bound = prev_ub[rxn.id]
    return growth, model.solver.status

def _gene_deletion(model, ids):
    all_reactions = []
    for g_id in ids:
        all_reactions.extend(
            find_gene_knockout_reactions(
                model, (model.genes.get_by_id(g_id),)
            )
        )
    growth, status = _reactions_knockouts_with_restore(model, all_reactions)
    return (ids, growth, status)

def extract_knockout_results(result_iter):
    result = pd.DataFrame([
        (frozenset(ids), growth, status)
        for (ids, growth, status) in result_iter
    ], columns=['ids', 'growth', 'status'])
    result.set_index('ids', inplace=True)
    return result

def _entities_ids(entities):
    try:
        return [e.id for e in entities]
    except AttributeError:
        return list(entities)

def _element_lists(entities, *ids):
    lists = list(ids)
    if lists[0] is None:
        lists[0] = entities
    result = [_entities_ids(lists[0])]
    for l in lists[1:]:
        if l is None:
            result.append(result[-1])
        else:
            result.append(_entities_ids(l))
    return result


def single_gene_knockout(model, gene_list = []):
    '''
    Context-specific programming causes the problem with TFA model.

    Parameters
    ----------
    model : cobra or pytfa model
        DESCRIPTION.
    gene_list : list of gene ids, optional
        DESCRIPTION. The default is [].

    Returns
    -------
    del_results : pd.DataFrame
        DESCRIPTION.

    '''
    if len(gene_list) == 0:
        gene_list = [x.id for x in model.genes]
    element_list = _element_lists(model.genes, gene_list)
    args = set([frozenset(comb) for comb in product(*element_list)])
    del_results = extract_knockout_results(map(
                partial(_gene_deletion, model), args))
        
    
    return del_results