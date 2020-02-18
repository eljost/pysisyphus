from pysisyphus.line_searches.LineSearch import LineSearch

class Dummy(LineSearch):

    def run_line_search(self):
        return 1.
