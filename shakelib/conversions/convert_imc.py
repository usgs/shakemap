from openquake.hazardlib import const


class ComponentConverter(object):
    """Base class for implementing conversions between components."""
    def __init__(self, imc_in, imc_out):
        self.imc_in = imc_in
        self.imc_out = imc_out
        graph = {'Average Horizontal (RotD50)': set([
                'Average Horizontal (GMRotI50)',
                'Average horizontal',
                'Horizontal Maximum Direction (RotD100)',
                'Greater of two horizontal']),
         'Average Horizontal (GMRotI50)': set([
                'Average Horizontal (RotD50)',
                'Greater of two horizontal']),
         'Average horizontal': set([
                'Average Horizontal (RotD50)',
                'Greater of two horizontal']),
         'Horizontal Maximum Direction (RotD100)': set([
                'Average Horizontal (RotD50)',
                'Greater of two horizontal']),
         'Greater of two horizontal': set([
                'Average Horizontal (RotD50)',
                'Average Horizontal (GMRotI50)',
                'Average horizontal',
                'Horizontal Maximum Direction (RotD100)'])}
        self.path = self.getShortestPath(graph, imc_in, imc_out)

    def getShortestPath(self, graph, start, goal):
        try:
            return next(self.bfs_paths(graph, start, goal))
        except StopIteration:
            return None

    def bfs_paths(self, graph, start, goal):
        queue = [(start, [start])]
        while queue:
            (vertex, path) = queue.pop(0)
            for next in graph[vertex] - set(path):
                if next == goal:
                    yield path + [next]
                else:
                    queue.append((next, path + [next]))
