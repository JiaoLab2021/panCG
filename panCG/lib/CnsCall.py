import re


def get_longest_trans_from_gff(gffFile, outFile):
    GeneDir = {}  # gene hash
    mRNADir = {}  # mRNA hash
    with open(gffFile, "r") as fi:
        # Get all genes
        while True:
            line = fi.readline()
            if not line:
                break
            if line.startswith("#") or not line.strip():  # Skip empty lines and comment lines
                continue
            line_li = line.strip().split("\t")
            if line_li[2] == "gene":
                geneID = re.findall(r'ID=([^;$]+)', line_li[8])[0]
                GeneDir[geneID] = {}
            elif line_li[2] == "mRNA":
                mRNAID = re.findall(r"ID=([^;$]+)", line_li[8])[0]
                mRNA_Parent = re.findall(r'Parent=([^;$]+)', line_li[8])[0]
                mRNADir[mRNAID] = line
                GeneDir[mRNA_Parent][mRNAID] = int(line_li[4]) - int(line_li[3])
            elif line_li[2].lower() == 'cds' or line_li[2].lower() == 'exon' or "utr" in line_li[2].lower():
                if len(re.findall(r'Parent=([^;$]+)', line_li[8])) == 0:
                    raise Exception("Missing Parent keyword")
                Parent = re.findall(r'Parent=([^;$]+)', line_li[8])[0]
                if Parent in mRNADir:
                    mRNADir[mRNAID] += line
                else:
                    raise Exception(f"{mRNAID} Missing Parent keyword")
            else:
                continue

    fo = open(outFile, "w")
    for gene_id, mRNA_length_dict in GeneDir.items():
        length = 0
        mRNAid = ""
        for mRNA_id, mRNA_length in mRNA_length_dict.items():
            if mRNA_length > length:
                length = mRNA_length
                mRNAid = mRNA_id
        fo.write(mRNADir[mRNAid])
    fo.close()
            

class SlideWindow:
    def __init__(self, chrID, size, high_PhastCons_score, low_PhastCons_score, trackDir, trackStart, trackEnd):
        self.chrID = chrID
        self.size = size
        self.high_PhastCons_score = high_PhastCons_score
        self.low_PhastCons_score = low_PhastCons_score
        self.start = trackStart
        self.end = trackStart + size - 1
        self.trackDir = trackDir
        self.trackStart = trackStart
        self.trackEnd = trackEnd
        self.trackDir[trackStart - 1] = 0
        self.trackDir[trackEnd + 1] = 0

    def judge(self):
        li = [self.trackDir[i] for i in range(self.start, self.end + 1)]
        if min(li) >= self.high_PhastCons_score:
            return True
        else:
            return False

    def extend_up(self):
        self.start = self.start - 1

    def extend_down(self):
        self.end = self.end + 1

    def extend_stop(self, bedFile):
        # print("{}\t{}\t{}\n".format(self.chrID, self.start, self.end))
        self.write_to_bed(bedFile)

    def WindowSlide(self):
        self.start = self.start + 1
        self.end = self.end + 1

    def NewWindow(self):
        self.start = self.end + 1
        self.end = self.start + self.size - 1

    def write_to_bed(self, bedFile):
        fi = open(bedFile, "a+")
        fi.write("{}\t{}\t{}\n".format(self.chrID, self.start - 1, self.end))  # Convert to bed format
        fi.close()

    @staticmethod
    def generateTrack(widFile):
        track_dict = {}
        start = 0
        n = 0
        chrID = ""
        fi = open(widFile, "r")
        while True:
            line = fi.readline()
            if not line:
                break
            line = line.strip()
            if "fixedStep" in line:
                chrID = re.findall(r"chrom=(.+?) ", line)[0]
                start = re.findall(r"start=(.+?) ", line)[0]
                step = re.findall(r"step=(.+?)", line)[0]
                if int(step) != 1:
                    raise ValueError("step != 1")
                n = int(start)
            else:
                track_dict[n] = float(line)
                n += 1
        fi.close()
        return track_dict, int(start), n - 1, chrID


def run_SlideWindow(high_PhastCons_score, low_PhastCons_score, wig, window, bedFile):
    track_dict, trackStart, trackEnd, chrID = SlideWindow.generateTrack(wig)
    SlideWindower = SlideWindow(chrID, window, high_PhastCons_score, low_PhastCons_score,
                                track_dict, trackStart, trackEnd)
    while True:
        if SlideWindower.end > SlideWindower.trackEnd:
            break
        if SlideWindower.judge():
            # Extend downward
            while True:
                if SlideWindower.trackDir[SlideWindower.end + 1] > SlideWindower.low_PhastCons_score:
                    # print(w.end)
                    SlideWindower.extend_down()
                else:
                    break
            # Extend upward
            while True:
                if SlideWindower.trackDir[SlideWindower.start - 1] > SlideWindower.low_PhastCons_score:
                    SlideWindower.extend_up()
                else:
                    break
            SlideWindower.extend_stop(bedFile)
            SlideWindower.NewWindow()
        else:
            SlideWindower.WindowSlide()

