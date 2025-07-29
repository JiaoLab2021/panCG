import pandas as pd
import re


class SyntenyAnchor:
    flank_length = 1000

    def __init__(self, simple_file):
        self.simple_file = simple_file

    def parse_simple(self):
        data = pd.read_csv(self.simple_file, sep="\t", header=None,
                           names=["que_anchor_1", "que_anchor_2", "ref_anchor_1", "ref_anchor_2", "score", "strand"])
        data[["que_chrID", "que_start", "que_end", "ref_chrID", "ref_start", "ref_end"]] = data.apply(self.add_block, axis=1)
        out_data = data[["que_chrID", "que_start", "que_end", "ref_chrID", "ref_start", "ref_end"]]
        out_data.loc[:, "que_start"] = out_data["que_start"] - 1000
        out_data.loc[:, "que_end"] = out_data["que_end"] + 1000
        out_data.loc[:, "ref_start"] = out_data["ref_start"] - 1000
        out_data.loc[:, "ref_end"] = out_data["ref_end"] + 1000
        print(out_data)

    @staticmethod
    def match_CNS_coordinate(CNSname):
        pattern = r'^(.*?):(\d+)-(\d+)'
        match = re.match(pattern, CNSname)
        if match:
            chrID, start, end = match.groups()
            return chrID, int(start), int(end)
        else:
            raise Exception(
                f"The format of {CNSname} is incorrect. It must be chrID:start-end. ':' cannot appear in chrID")

    def add_block(self, row):
        que_chrID_1, que_start_1, que_end_1 = self.match_CNS_coordinate(row["que_anchor_1"])
        que_chrID_2, que_start_2, que_end_2 = self.match_CNS_coordinate(row["que_anchor_2"])
        if que_chrID_1 == que_chrID_2:
            que_chrID = que_chrID_1
            que_start = min(que_start_1, que_end_1, que_start_2, que_end_2)
            que_end = max(que_start_1, que_end_1, que_start_2, que_end_2)
        else:
            raise Exception("The chrID number is wrong")

        ref_chrID_1, ref_start_1, ref_end_1 = self.match_CNS_coordinate(row["ref_anchor_1"])
        ref_chrID_2, ref_start_2, ref_end_2 = self.match_CNS_coordinate(row["ref_anchor_2"])
        if ref_chrID_1 == ref_chrID_2:
            ref_chrID = ref_chrID_1
            ref_start = min(ref_start_1, ref_end_1, ref_start_2, ref_end_2)
            ref_end = max(ref_start_1, ref_end_1, ref_start_2, ref_end_2)
        else:
            raise Exception("The chrID number is wrong")

        return pd.Series([que_chrID, que_start, que_end, ref_chrID, ref_start, ref_end])


pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_columns', None)

simple_file = "/Users/ltan/Documents/tmp/3/test.simple"
SyntenyAnchorer = SyntenyAnchor(simple_file)
SyntenyAnchorer.parse_simple()



