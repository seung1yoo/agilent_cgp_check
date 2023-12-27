

import sys
from pathlib import Path

class ExtractGeneInfoFromPanelBed:
    def __init__(self):
        self.bg_dic = dict()
        self.panel_dic = dict()

    def make_background(self, file_path):

        for line in file_path.open():
            items = line.rstrip("\n").split("\t")
            if items[0] in ["Gene stable ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            gene_id = items[idx_dic["Gene name"]]
            ens_id = items[idx_dic["Gene stable ID"]]
            if gene_id:
                self.bg_dic.setdefault(gene_id, {})
                self.bg_dic[gene_id].setdefault("ens_id", ens_id)
            else:
                pass

        return 1

    def update_panel_info(self, file_path):

        for line in file_path.open():
            items = line.rstrip("\n").split("\t")
            if len(items) not in [4]:
                print(items)
                continue
            _gene_id = items[3]
            if ";" in _gene_id:
                for gene_id in _gene_id.split(";"):
                    self.panel_dic.setdefault(gene_id, {})
                    self.add_count_of_probes(gene_id)
                    self.add_covered_position(gene_id, items)
            elif "," in _gene_id:
                for gene_id in _gene_id.split(","):
                    self.panel_dic.setdefault(gene_id, {})
                    self.add_count_of_probes(gene_id)
                    self.add_covered_position(gene_id, items)
            else:
                self.panel_dic.setdefault(_gene_id, {})
                self.add_count_of_probes(_gene_id)
                self.add_covered_position(_gene_id, items)

    def add_count_of_probes(self, gene_id):
        self.panel_dic[gene_id].setdefault("n.probes", 0)
        self.panel_dic[gene_id]["n.probes"] += 1

    def add_covered_position(self, gene_id, items):
        self.panel_dic[gene_id].setdefault("covered.position", {})
        for position in range(int(items[1]), int(items[2])+1):
            _key = f"{items[0]}:{position}"
            self.panel_dic[gene_id]["covered.position"].setdefault(_key, 0)
            self.panel_dic[gene_id]["covered.position"][_key] += 1

    def write_table(self, outfn):
        outfh = Path(outfn).open("w")
        headers = ["gene_id"]
        headers.append("ens_id")
        headers.append("n.probes")
        headers.append("bp.covered.region.by.probes.1")
        headers.append("bp.covered.region.by.probes.2")
        outfh.write("{0}\n".format("\t".join(headers)))
        for gene_id, info_dic in sorted(self.panel_dic.items()):
            items = [gene_id]
            if gene_id in self.bg_dic:
                items.append(self.bg_dic[gene_id]["ens_id"])
            else:
                items.append("-")
            items.append(str(info_dic["n.probes"]))
            #
            items.append(str(self.cal_covered_length(info_dic["covered.position"], 1)))
            items.append(str(self.cal_covered_length(info_dic["covered.position"], 2)))
            outfh.write("{0}\n".format("\t".join(items)))
        outfh.close()

    def cal_covered_length(self, _dic, min_cov_depth):
        count = 0
        for _key, cov_depth in _dic.items():
            if cov_depth >= min_cov_depth:
                count += 1
        return count


    def print_target_gene_count(self):
        print(sorted([x for x in self.panel_dic]))
        print(len([x for x in self.panel_dic if not x in ["CNVBackbone", "MSICovered"]]))

    def print_panel_size(self):
        panel_size = 0
        for gene_id, info_dic in sorted(self.panel_dic.items()):
            panel_size += self.cal_covered_length(info_dic["covered.position"], 1)
        print(panel_size)

def main():
    obj = ExtractGeneInfoFromPanelBed()
    obj.make_background(Path("hg38_biomart_ens_id_map.txt"))
    obj.update_panel_info(Path("../A3416642_Covered.bed"))
    obj.write_table("cgp_gene_info.tsv")
    obj.print_target_gene_count()
    obj.print_panel_size()




if __name__=="__main__":
    main()
