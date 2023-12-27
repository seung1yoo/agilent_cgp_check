
from pathlib import Path
import sys


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

def extract_anno_dic(anno):
    anno_dic = dict()
    for _anno in anno.split(";"):
        try:
            key, value = _anno.split(' "')
            #print(key.strip(), value.strip('"'))
            anno_dic.setdefault(key.strip(), value.strip('"'))
        except ValueError as e:
            print(f"def extract_gene_id : {_anno} {e}")
            continue
    return anno_dic



def main():
    obj = ExtractGeneInfoFromPanelBed()
    obj.make_background(Path("extract_target_gene_list/hg38_biomart_ens_id_map.txt"))
    obj.update_panel_info(Path("A3416642_Covered.bed"))



    outfn = "MANE.GRCh38.v1.2.refseq_genomic.exon.CoverageWith_A3416642.tsv"
    outfh = open(outfn, "w")
    headers = ["chr"]
    headers.append("start")
    headers.append("end")
    headers.append("strand")
    headers.append("gene_id")
    headers.append("ens_tid")
    headers.append("nm_id")
    headers.append("n.probes")
    headers.append("covered_region_by_probes")
    headers.append("exonic_region")
    headers.append("ratio_of_covered")
    headers.append("total.probes_of_this_gene")
    headers.append("bases.covered_region_by_all_probes_on_this_gene")
    outfh.write("{0}\n".format("\t".join(headers)))
    for line in open("MANE.GRCh38.v1.2.refseq_genomic.exon.CoverageWith_A3416642.bed"):
        items = line.rstrip("\n").split("\t")
        new_items = [items[0]] #chr
        new_items.append(items[3]) #start
        new_items.append(items[4]) #end
        new_items.append(items[6]) #strand
        new_items.append(extract_anno_dic(items[8])["gene_id"]) #gene_id
        new_items.append(extract_anno_dic(items[8])["db_xref"]) #ens_tid
        new_items.append(extract_anno_dic(items[8])["transcript_id"]) #nm_id
        new_items.append(items[9]) #n.probe
        new_items.append(items[10]) #covered_region_by_probe
        new_items.append(items[11]) #exonic_region
        new_items.append(items[12]) #ratio_of_covered

        gene_id = extract_anno_dic(items[8])["gene_id"]
        if gene_id in obj.panel_dic:
            new_items.append(obj.panel_dic[gene_id]["n.probes"]) #total.probes_of_this_gene
            new_items.append(obj.cal_covered_length(obj.panel_dic[gene_id]["covered.position"], 1)) #total.covered.region_by_all_probes_of_this_gene
        else:
            new_items.append("-")
            new_items.append("-")
        outfh.write("{0}\n".format("\t".join([str(x) for x in new_items])))
    outfh.close()


if __name__=='__main__':
    main()
