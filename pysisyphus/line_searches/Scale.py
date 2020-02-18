from pysisyphus.line_searches.LineSearch import LineSearch

class Scale(LineSearch):

    def run_line_search(self):
        return 1.
