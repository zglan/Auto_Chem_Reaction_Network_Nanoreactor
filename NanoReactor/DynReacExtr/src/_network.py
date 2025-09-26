from collections import Counter

import cmocean
import networkx as nx
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

""" continuous colorbar  """
# cmap = cm.Miscellaneous.reversed()
# cmap = cm.plasma.reversed()

""" discrete intervals colorbar """
# cmap = (mpl.colors.ListedColormap(
#     ['royalblue', 'cyan', 'yellow', 'orange'])
#          .with_extremes(over='0.25', under='0.75'))

""" pos setting """
# pos = nx.spring_layout(G, k=0.5, iterations=200, scale=0.5)
# pos = nx.kamada_kawai_layout(G)

""" norm setting """
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='max')


default_net_param = {
    'node_size': 600,
    'r_node_size': 3,

    'edge_width': 2,
    'r_edge_width': 100,

    'r_alpha': 12,
    'node_alpha': 0.5,
    'edge_alpha': 0.15,
    'node_alpha_criter': 0.6,
    'edge_alpha_criter': 0.7,

    'bound': [3, 6, 10, 20, 35, 50],
    
    'cbar_label_size': 20,
    'cbar_tick_size': 14,
    'cbar_label_font': 'Detetected Times',

    'net_label_size': 16,
    'spec_label': {1: '1'},
    'spec_label_color': 'Gray',

    'fig_size': [24, 15],
    'fig_dpi': 400,
    'fig_name': 'reac_Event_Network.png'
    }

def upadteParam(default_param=default_net_param, new_param={}):
    default_param.update(new_param)
    return default_param

class ReacEveNet:

    def __init__(self, smi_dict, smi_list, idx_rela_list, net_param=default_net_param):
        self.node_list = [smi_dict[smi] for smi in Counter(smi_list).keys()]
        self.edge_list = list(Counter(idx_rela_list).keys())

        self.node_wt = list(Counter(smi_list).values())
        self.edge_wt = list(Counter(idx_rela_list).values())
        self.tot_wt = len(idx_rela_list)

        for key in net_param:
            setattr(self, key, net_param[key])
        
        self.norm = mcolors.Normalize(vmin=3, vmax=50)
        self.cmap = cmocean.cm.haline
        self.layout_param = {}
        self.layout_method = nx.spring_layout

    def __setAlpha(self, wt, base, criter):
        alpha = wt / self.tot_wt * self.r_alpha + base
        if alpha <= criter:
            return alpha
        else:
            return 1
    
    def __setSize(self, wt, base, r):
        return wt * r + base

    def __setAttr(self, wt_list, 
            size, r_size, 
            alpha, alpha_criter):
        color_list = wt_list
        size_list = [self.__setSize(wt, size, r_size)
            for wt in wt_list]
        alpha_list = [self.__setAlpha(wt, alpha, alpha_criter)
            for wt in wt_list]
        return color_list, size_list, alpha_list
    
    def layout_method(self, *argc):
        self.layout_method(*argc)
    
    def __buildDiNet(self):
        G = nx.DiGraph()
        G.add_nodes_from(self.node_list)
        G.add_edges_from(self.edge_list)
        return G

    def __setPos(self, G):
        params = {}
        params['G'] = G
        for key in self.layout_param.keys():
            params[key] = self.layout_param[key]
        pos = self.layout_method(**params)
        return pos
        
    def __setNodeEdge(self, G, pos):
        node_color_list, node_size_list, \
            node_alpha_list = self.__setAttr(self.node_wt, 
                self.node_size, self.r_node_size,
                self.node_alpha, self.node_alpha_criter)

        edge_color_list, edge_width_list, \
            edge_alpha_list = self.__setAttr(self.edge_wt, 
                self.edge_width, self.r_edge_width / self.tot_wt,
                self.edge_alpha, self.edge_alpha_criter)

        nodes = nx.draw_networkx_nodes(
            G, pos, 
            nodelist=self.node_list, 
            node_size=node_size_list, 
            node_color=node_color_list,
            alpha=node_alpha_list,
            cmap=self.cmap)
        edges = nx.draw_networkx_edges(
            G, pos,
            edgelist=self.edge_list,
            width=edge_width_list,
            edge_color=edge_color_list, 
            alpha=edge_alpha_list,
            edge_cmap=self.cmap,
            arrowstyle="->", 
            arrowsize=25,
            connectionstyle="arc3,rad=0.4")
        labels = nx.draw_networkx_labels(
            G, pos, 
            font_size=self.net_label_size, 
            font_family="sans-serif", 
            font_color='w',
            clip_on=False)
        spec_labels = nx.draw_networkx_labels(
            G, pos, 
            self.spec_label, 
            font_size=self.net_label_size, 
            font_family="sans-serif",
            font_color=self.spec_label_color)

        return G
    
    def __setColorBar(self, ax):
        color_bar = plt.colorbar(
            cm.ScalarMappable(
                cmap=self.cmap, norm=self.norm),
            ticks=self.bound,
            extend='max',
            ax=ax,
            pad=0)
            # orientation='horizontal',
            # spacing='proportional',
        color_bar.ax.tick_params(
            labelsize=self.cbar_tick_size)
        color_bar.set_label(
            self.cbar_label_font, 
            size=self.cbar_label_size,
            loc='center')
        
        return color_bar
    
    def setLayoutParam(self, param={}):
        self.layout_param = param

    def setLayoutMethod(self, lm):
        self.layout_method = lm

    def setCmap(self, cmap):
        self.cmap = cmap

    def setNorm(self, norm):
        self.norm = norm

    def drawNet(self):
        G = self.__buildDiNet()
        pos = self.__setPos(G)
        G = self.__setNodeEdge(G, pos)

        ax = plt.gca()
        ax.set_axis_off()
        cb = self.__setColorBar(ax)

        fig = plt.gcf()
        fig.set_size_inches(
            self.fig_size[0], self.fig_size[1])
        
        plt.savefig(
            self.fig_name, 
            dpi=self.fig_dpi, 
            bbox_inches='tight')

