def safe_mean(arr):
    """安全均值计算，出错返回 None"""
    try:
        if arr is None or len(arr) == 0:
            return None
        return float(np.nanmean(arr))
    except Exception:
        return None


def parse_confidences_json(conf_path, pdb_path):
    try:
        with open(conf_path) as f:
            conf = json.load(f)

        chains = get_chains_from_pdb(pdb_path)
        token_chain_ids, token_res_ids = extract_token_chain_and_res_ids(pdb_path)
        pae = np.array(conf.get("pae", []))

        # 长度一致性检查
        if len(token_chain_ids) != pae.shape[0]:
            return None, None, None

        # 计算每条链的 PAE
        chain_indices = {chain: [] for chain in chains}
        for i, chain in enumerate(token_chain_ids):
            chain_indices[chain].append(i)

        chain_pae = {
            chain: safe_mean(pae[np.ix_(idxs, idxs)]) if idxs else None
            for chain, idxs in chain_indices.items()
        }

        ipae = {}
        pae_interaction = {}

        for ch1, ch2 in combinations(chains, 2):
            try:
                idx1_, idx2_ = get_interface_res_from_pdb(pdb_path, chain1=ch1, chain2=ch2)
                if not idx1_ or not idx2_:
                    continue

                idx1 = [i for i, (res_id, chain) in enumerate(zip(token_res_ids, token_chain_ids))
                        if chain == ch1 and res_id in idx1_]
                idx2 = [i for i, (res_id, chain) in enumerate(zip(token_res_ids, token_chain_ids))
                        if chain == ch2 and res_id in idx2_]

                if idx1 and idx2:
                    ipae[f"{ch1}_{ch2}"] = safe_mean([
                        safe_mean(pae[np.ix_(idx1, idx2)]),
                        safe_mean(pae[np.ix_(idx2, idx1)])
                    ])

                chain_1_indices = [i for i, c in enumerate(token_chain_ids) if c == ch1]
                chain_2_indices = [i for i, c in enumerate(token_chain_ids) if c == ch2]

                pae_interaction[f"{ch1}_{ch2}"] = safe_mean([
                    safe_mean(pae[np.ix_(chain_1_indices, chain_2_indices)]),
                    safe_mean(pae[np.ix_(chain_2_indices, chain_1_indices)])
                ])

            except Exception as e:
                print(f"[Warning] Failed to process pair ({ch1}, {ch2}): {e}")
                ipae[f"{ch1}_{ch2}"] = None
                pae_interaction[f"{ch1}_{ch2}"] = None

        return chain_pae, ipae, pae_interaction
    except Exception as e:
        print(f"[Warning] Failed in parse_confidences_json: {e}")
        return None, None, None


def process_single_description(args):
    description, input_pdb_dir, base_dir = args
    try:
        base_path = Path(base_dir) / description / "seed-10_sample-0"
        summary_path = base_path / "summary_confidences.json"
        conf_path = base_path / "confidences.json"
        pdb_path = Path(input_pdb_dir) / f"{description}.pdb"

        if not (summary_path.exists() and pdb_path.exists() and conf_path.exists()):
            return None, f"{description}: missing files"

        summary = json.loads(summary_path.read_text())
        conf = json.loads(conf_path.read_text())
        chains = get_chains_from_pdb(pdb_path)

        iptm = dict(zip(chains, summary.get("chain_iptm", [])))
        ptm = dict(zip(chains, summary.get("chain_ptm", [])))

        interchain_iptm_dict = {}
        iptm_matrix = summary.get("chain_pair_iptm", [])
        num_chains = len(chains)
        for i in range(num_chains):
            for j in range(i + 1, num_chains):
                try:
                    interchain_iptm_dict[f"{chains[i]}_{chains[j]}"] = iptm_matrix[i][j]
                except Exception:
                    interchain_iptm_dict[f"{chains[i]}_{chains[j]}"] = None

        atom_plddts = conf.get("atom_plddts", [])
        atom_chain_ids = conf.get("atom_chain_ids", [])
        chain_plddt = {}
        for ch in chains:
            chain_plddt[ch] = safe_mean([pl for pl, cid in zip(atom_plddts, atom_chain_ids) if cid == ch])
        complex_plddt = safe_mean(list(chain_plddt.values()))

        try:
            chain_pae, ipae, inter_pae = parse_confidences_json(conf_path, str(pdb_path))
        except Exception:
            chain_pae, ipae, inter_pae = None, None, None

        return {
            "description": description,
            "AF3Score_chain_ptm": ptm,
            "AF3Score_chain_iptm": iptm,
            "AF3Score_interchain_iptm": interchain_iptm_dict,
            "AF3Score_complex_plddt": complex_plddt,
            "AF3Score_chain_ca_plddt": chain_plddt,
            "AF3Score_chain_pae": chain_pae,
            "AF3Score_ipae": ipae,
            "AF3Score_pae_interaction": inter_pae,
        }, None

    except Exception as e:
        return None, f"{description}: {str(e)}"
