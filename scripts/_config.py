# -*- coding: utf-8 -*-
"""Single configuration loader for every Python script in this repository.

Paths, the sample list and the analysis thresholds all come from config/config.ini,
so that no script contains a machine-specific path. Import it as:

    from _config import cfg
    cfg.results_dir          # absolute path to results/
    cfg.samples              # ['MOCK_1', ..., 'PR8_3']
    cfg.is_mock              # [True, True, True, False, False, False]
    cfg.get('analysis', 'fdr_cut', 0.05)
"""
import configparser
import os

# repository root = the parent of the directory holding this file (scripts/)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_PATH = os.path.join(ROOT, "config", "config.ini")


def _split(value):
    return [v.strip() for v in value.split(",") if v.strip()]


class Config:
    def __init__(self, path=CONFIG_PATH):
        if not os.path.exists(path):
            raise SystemExit("STOP: configuration file not found: %s" % path)
        self.parser = configparser.ConfigParser(inline_comment_prefixes=("#",))
        self.parser.read(path, encoding="utf-8")
        p, s = self.parser["paths"], self.parser["samples"]

        # the repository root, so that scripts can reach docs/ and other non-data paths
        self.root = ROOT

        # data_root may be absolute or relative to the repository root
        self.data_root = p.get("data_root", "data")
        if not os.path.isabs(self.data_root):
            self.data_root = os.path.join(ROOT, self.data_root)
        self.data_root = os.path.normpath(self.data_root)

        def under(root, key, default):
            v = p.get(key, default)
            return v if os.path.isabs(v) else os.path.normpath(os.path.join(root, v))

        self.raw_dir = under(self.data_root, "raw_dir", "raw")
        self.processed_dir = under(self.data_root, "processed_dir", "processed")
        self.qc_dir = under(self.data_root, "qc_dir", "qc_evidence")
        self.interim_dir = under(self.data_root, "interim_dir", "interim")
        self.results_dir = under(ROOT, "results_dir", "results")
        self.figdata_dir = under(ROOT, "figdata_dir", "figdata")
        self.figures_dir = under(ROOT, "figures_dir", "figures")
        self.tables_dir = under(ROOT, "tables_dir", "tables")

        self.samples = _split(s.get("names", ""))
        self.groups = _split(s.get("groups", ""))
        self.reference = s.get("reference", "MOCK")
        self.case = s.get("case", "PR8")
        if len(self.samples) != len(self.groups):
            raise SystemExit("STOP: [samples] names and groups differ in length")
        self.is_mock = [g == self.reference for g in self.groups]
        self.is_case = [g == self.case for g in self.groups]

        self.geo = dict(self.parser["geo"])
        self.gsm = dict(zip(self.samples, _split(self.parser["geo"].get("gsm", ""))))
        self.input_files = dict(self.parser["input_files"])
        self.fe_files = _split(self.input_files.get("fe_files", ""))
        self.qc_files = dict(self.parser["qc_files"]) if self.parser.has_section("qc_files") else {}

    # -- helpers ---------------------------------------------------------
    def get(self, section, option, default=None):
        if not self.parser.has_option(section, option):
            return default
        return self.parser.get(section, option)

    def getfloat(self, section, option, default):
        return self.parser.getfloat(section, option) if self.parser.has_option(section, option) else default

    def getint(self, section, option, default):
        return self.parser.getint(section, option) if self.parser.has_option(section, option) else default

    def raw(self, name):
        return os.path.join(self.raw_dir, name)

    def deposit(self, name):
        """Locate a deposited file by name: raw/ first, then processed/.

        Matrix.xlsx lives in raw/ with the Feature Extraction files; the provider
        result tables live in processed/. Callers should not have to know which.
        """
        for d in (self.raw_dir, self.processed_dir):
            p = os.path.join(d, name)
            if os.path.exists(p):
                return p
        return os.path.join(self.raw_dir, name)

    def processed(self, key):
        return os.path.join(self.processed_dir, self.input_files[key])

    def qc(self, key):
        """Provider QC report shipped with the repository (not part of the GEO deposit)."""
        if key not in self.qc_files:
            raise SystemExit("STOP: no file named for qc_files.%s in config/config.ini" % key)
        return os.path.join(self.qc_dir, self.qc_files[key])

    def interim(self, name, create=False):
        return self._mk(self.interim_dir, name, create)

    def result(self, name, create=False):
        return self._mk(self.results_dir, name, create)

    def figdata(self, name, create=False):
        return self._mk(self.figdata_dir, name, create)

    def figure(self, name, create=False):
        return self._mk(self.figures_dir, name, create)

    def table(self, name, create=False):
        return self._mk(self.tables_dir, name, create)

    @staticmethod
    def _mk(d, name, create):
        if create:
            os.makedirs(d, exist_ok=True)
        return os.path.join(d, name)

    def make_dirs(self):
        for d in (self.interim_dir, self.results_dir, self.figdata_dir, self.figures_dir, self.tables_dir):
            os.makedirs(d, exist_ok=True)

    def check_inputs(self):
        """Report any expected input file that is missing, so a run fails early and clearly."""
        want = [self.raw(self.input_files["matrix"])] + [self.raw(n) for n in self.fe_files]
        want += [self.processed(k) for k in ("mrna_q_flagpass", "mrna_expr_flagpass",
                                             "lncrna_q_flagpass", "others_q_flagpass")]
        missing = [w for w in want if not os.path.exists(w)]
        if missing:
            print("missing input files (%d):" % len(missing))
            for m in missing:
                print("  -", m)
        return missing

    def report(self):
        print("repository root : %s" % ROOT)
        print("data root       : %s" % self.data_root)
        print("raw / processed : %s | %s" % (self.raw_dir, self.processed_dir))
        print("samples         : %s" % ", ".join("%s(%s)" % (s, g) for s, g in zip(self.samples, self.groups)))
        print("thresholds      : FDR<%g | log2FC>=%g | Dpct>=%g pp" % (
            self.getfloat("analysis", "fdr_cut", 0.05),
            self.getfloat("analysis", "log2fc_cut", 1.0),
            self.getfloat("analysis", "pctmod_pp_cut", 5.0)))


cfg = Config()

if __name__ == "__main__":
    cfg.report()
    missing = [p for p in ([cfg.raw(n) for n in ["Matrix.xlsx"] + cfg.fe_files] +
                           [cfg.processed(k) for k in
                            ("mrna_q_flagpass", "mrna_expr_flagpass", "lncrna_q_flagpass", "others_q_flagpass")])
               if not os.path.exists(p)]
    print("\nmissing input files: %d" % len(missing))
    for m in missing:
        print("  -", m)
