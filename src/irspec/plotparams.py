import matplotlib.pyplot as plt 
import matplotlib

class PlotParams:
    
    
    def __init__(self, palatte="dark", scaling="poster", font="Helvetica"):
        """ 
        Initializes an instance of the PlotParams object. When called,
        resets the global plotting parameters for any script.
        
        Arguments
        ---------
        style : str
            The broad classification for what style is to be used in 
            figures
        font : str
            The font family to be used in figures
        """
        palattes = ["light", "dark"]
        scalings = ["paper", "presentation", "poster"]
        if palatte not in palattes:
            raise ValueError("Unsupported palatte was implemented")
        self.palatte = palatte
        if scaling not in scalings:
            raise ValueError("Unsupported scaling was implemented")
        self.scaling = scaling
        
        if self.palatte == "dark":
            self.apply_dark()
        if self.palatte == "light":
            self.apply_light()
        
        if self.scaling == "poster":
            self.apply_poster_scaling()
        if self.scaling == "presentation":
            self.apply_presentation_scaling()
        if self.scaling == "paper":
            self.apply_paper_scaling()
        
        plt.rcParams["font.family"] = font
    
    
    def hide_axes_labels(self):
        plt.rcParams["axes.labelsize"] = 0
        plt.rcParams["xtick.labelsize"] = 0
        plt.rcParams["ytick.labelsize"] = 0
        plt.rcParams["axes.grid"] = False
    
    def apply_poster_scaling(self):
        
        plt.rcParams["axes.labelsize"] = 32
        
        plt.rcParams["xtick.labelsize"] = 0
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.labelsize"] = 0
        plt.rcParams["ytick.direction"] = "in"
        
        plt.rcParams["axes.grid"] = True
        plt.rcParams["grid.alpha"] = 0.5
        plt.rcParams["grid.linewidth"] = 0.5
        plt.rcParams["legend.fontsize"] = 12
        
        plt.rcParams["savefig.dpi"] = 1200
        plt.rcParams["savefig.bbox"] = "tight"
    
    
    def apply_presentation_scaling(self):
        
        plt.rcParams["axes.labelsize"] = 22
        
        plt.rcParams["xtick.labelsize"] = 20
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.labelsize"] = 20
        plt.rcParams["ytick.direction"] = "in"
        
        plt.rcParams["axes.grid"] = True
        plt.rcParams["grid.alpha"] = 0.5
        plt.rcParams["grid.linewidth"] = 0.5
        plt.rcParams["legend.fontsize"] = 12
        
        plt.rcParams["savefig.dpi"] = 1200
        plt.rcParams["savefig.bbox"] = "tight"
    
    
    def apply_paper_scaling(self):
        
        plt.rcParams["axes.labelsize"] = 16
        
        plt.rcParams["xtick.labelsize"] = 14
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.labelsize"] = 14
        plt.rcParams["ytick.direction"] = "in"
        
        #plt.rcParams["axes.grid"] = True
        #plt.rcParams["grid.alpha"] = 0.5
        #plt.rcParams["grid.linewidth"] = 0.5
        plt.rcParams["legend.fontsize"] = 12
        
        plt.rcParams["savefig.dpi"] = 1200
        plt.rcParams["savefig.bbox"] = "tight"
    
    
    def apply_light(self):
        """ 
        Applies dark mode settings with publication quality font.
        """
        
        plt.rcParams["figure.figsize"] = (8,6)
        
        
        plt.rcParams["axes.labelsize"] = 32
        plt.rcParams["axes.titlesize"] = 28
        
        plt.rcParams["xtick.labelsize"] = 28
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.labelsize"] = 28
        plt.rcParams["ytick.direction"] = "in"
        
        plt.rcParams["axes.grid"] = True
        plt.rcParams["grid.alpha"] = 0.5
        plt.rcParams["grid.linewidth"] = 0.5
        plt.rcParams["legend.fontsize"] = 12
        
        plt.rcParams["savefig.dpi"] = 800
        plt.rcParams["savefig.bbox"] = "tight"
    
    
    def apply_dark(self):
        """ 
        Applies dark mode settings with publication quality font.
        """
        
        plt.style.use('dark_background')
        plt.rcParams["figure.figsize"] = (8,6)
        
        plt.rcParams["axes.labelsize"] = 32
        plt.rcParams["axes.titlesize"] = 28
        
        plt.rcParams["xtick.labelsize"] = 28
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.labelsize"] = 28
        plt.rcParams["ytick.direction"] = "in"
        
        #plt.rcParams["axes.grid"] = False
        #plt.rcParams["grid.alpha"] = 0.5
        plt.rcParams["grid.linewidth"] = 0.5
        plt.rcParams["legend.fontsize"] = 12
        
        plt.rcParams["savefig.dpi"] = 800
        plt.rcParams["savefig.bbox"] = "tight"
    
    
    def dark_colors(self):
        return ["#e60049", "#0bb4ff", "#ffa300", "#50e991", "#9b19f5", "#e6d800", "#dc0ab4", "#b3d4ff", "#00bfa0"]
    
    def light_colors(self):
        return plt.get_cmap('tab20')

    def linestyles(self):
        return list(matplotlib.lines.lineStyles.keys())

    def markers(self):
        return list(matplotlib.markers.MarkerStyle.markers.keys())
