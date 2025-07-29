"""
CNS identification strategy:
    Identification based on scores in wig files.
    Use <window> bp as the window, slide according to the genome track, and stop sliding when all windows with scores greater than <high_PhastCons_score> are encountered.
    At this time, the extension is extended to both sides, and the extension is terminated when the score is less than <low_PhastCons_score>.
"""

import argparse
import re


class Window:
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


def get_parser():
    parser = argparse.ArgumentParser(
        description="This script extracts the conserved sequence in chromosome according to the score in WIG file.")

    parser.add_argument("-w", "--wig",
                        dest="wig",
                        help="Input wig format file.",
                        required=True)
    parser.add_argument("-high", "--high_PhastCons_score",
                        dest="high_PhastCons_score",
                        help="The high PhastCons score.",
                        required=True,
                        type=float)
    parser.add_argument("-low", "--low_PhastCons_score",
                        dest="low_PhastCons_score",
                        help="The low PhastCons score.",
                        required=True,
                        type=float)
    parser.add_argument("-window", "--window",
                        dest="window",
                        help="The sliding window size",
                        required=True,
                        type=int)
    parser.add_argument("-o", "--out_bed_file",
                        dest="out_bed_file",
                        help="The output bed file.",
                        required=True)
    return parser


def generateTrack(widFile):
    trackDir = {}
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
            trackDir[n] = float(line)
            n += 1
    fi.close()
    return trackDir, int(start), n - 1, chrID


def main():
    # high_PhastCons_score = 0.82
    # low_PhastCons_score = 0.55
    # window = 8
    # wig = "1.wig"
    # trackDir, trackStart, trackEnd, chrID = generateTrack(wig)
    # bedFile = "test.bed"
    # w = Window(chrID, window, high_PhastCons_score, low_PhastCons_score, trackDir, trackStart, trackEnd)

    parser = get_parser()
    args = parser.parse_args()
    high_PhastCons_score = args.high_PhastCons_score
    low_PhastCons_score = args.low_PhastCons_score
    wig = args.wig
    window = args.window
    bedFile = args.out_bed_file
    trackDir, trackStart, trackEnd, chrID = generateTrack(wig)
    w = Window(chrID, window, high_PhastCons_score, low_PhastCons_score, trackDir, trackStart, trackEnd)

    while True:
        if w.end > w.trackEnd:
            break
        if w.judge():
            # Extend downward
            while True:
                if w.trackDir[w.end + 1] > w.low_PhastCons_score:
                    # print(w.end)
                    w.extend_down()
                else:
                    break
            # Extend upward
            while True:
                if w.trackDir[w.start - 1] > w.low_PhastCons_score:
                    w.extend_up()
                else:
                    break
            w.extend_stop(bedFile)
            w.NewWindow()
        else:
            w.WindowSlide()


if __name__ == '__main__':
    main()

