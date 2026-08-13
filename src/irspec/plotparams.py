import matplotlib.pyplot as plt
import matplotlib

class PlotParams:


    # The supported palettes, and the one modules should reach for when the
    # caller expresses no preference. Light is the default because the maps
    # are primarily made for papers; dark is the presentation/poster variant.
    PALATTES = ("light", "dark")
    SCALINGS = ("paper", "presentation", "poster")
    DEFAULT_PALATTE = "light"
    DEFAULT_SCALING = "poster"

    # The house font. Declared here rather than left as a literal in the
    # signature so it can be referred to, and changed, in one place.
    DEFAULT_FONT = "Helvetica"

    # What the font degrades to where the house font is not installed. Without
    # a list matplotlib falls back to DejaVu Sans, which is wider and rounder
    # and does not sit well beside the rest of these figures.
    FONT_FALLBACKS = ("Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans")

    # Math is *not* set in the house font. macOS ships Helvetica without Greek
    # coverage, so aiming mathtext straight at it (`mathtext.fontset =
    # "custom"`) silently renders \chi as a Latin "X" and \nu as a "v" -- no
    # warning, no missing-glyph box -- which quietly corrupts a label like
    # $\chi^2_\nu$ into "X^2_v". STIX Sans is a sans-serif math face built to
    # sit beside Helvetica-class text and it covers the Greek letters and
    # italic variables these figures use.
    MATHTEXT_FONTSET = "stixsans"

    # Palette-dependent drawing standards. Anything that has to contrast
    # against the figure background -- overlay annotations, the scalebar,
    # the color a colormap paints NaN/masked spaxels -- resolves through
    # here rather than hardcoding "black" or "white" at the call site, so a
    # figure stays legible in whichever palette it was rendered under.
    _STANDARDS = {
        "light": {"foreground": "black", "background": "white",
                  "bad": "white"},
        "dark":  {"foreground": "white", "background": "black",
                  "bad": "black"},
    }

    # Every text size a figure draws, in points, one row per scaling. Sizes
    # that used to be hardcoded at the call sites -- colorbar labels, the
    # scalebar, map annotations -- read `annotation` and `axes_label` from
    # here, so switching scaling moves all of them together instead of
    # leaving 24pt text beside 14pt axis labels.
    #
    # `annotation` sits a step below `axes_label`: text drawn inside the axes
    # is a caption on the data, and should not compete with the frame.
    _SCALING_FONTS = {
        "paper":        {"axes_label": 16, "title": 18, "tick": 14,
                         "annotation": 14, "legend": 12},
        "presentation": {"axes_label": 22, "title": 24, "tick": 20,
                         "annotation": 20, "legend": 16},
        "poster":       {"axes_label": 32, "title": 34, "tick": 28,
                         "annotation": 28, "legend": 22},
    }

    def __init__(self, palatte=DEFAULT_PALATTE, scaling=DEFAULT_SCALING,
                 font=DEFAULT_FONT, grid=False):
        """
        Initializes an instance of the PlotParams object. When called,
        resets the global plotting parameters for any script.

        Arguments
        ---------
        palatte : str
            Either ``light`` or ``dark``; selects the figure background and
            the contrasting foreground standards (default ``light``)
        scaling : str
            The broad classification for what style is to be used in
            figures
        font : str
            The font family to be used in figures (default ``Helvetica``;
            see ``FONT_FALLBACKS`` for what it degrades to, and
            ``MATHTEXT_FONTSET`` for why math is set separately)
        grid : bool
            Draw the axes grid (celestial gridlines on a WCS projection).
            Off by default, which suits the spaxel maps -- a grid over an
            image obscures the data. Scripts that plot curves rather than
            images pass ``grid=True``.
        """
        if palatte not in self.PALATTES:
            raise ValueError(
                f"Unsupported palatte {palatte!r}; expected one of "
                f"{list(self.PALATTES)}")
        self.palatte = palatte
        if scaling not in self.SCALINGS:
            raise ValueError(
                f"Unsupported scaling {scaling!r}; expected one of "
                f"{list(self.SCALINGS)}")
        self.scaling = scaling
        self.font = font
        self.grid = grid

        self.apply()


    def apply(self):
        """Asserts this instance's palette, scaling, and font onto the global
        rcParams.

        rcParams are process-global, so any other module that styles a figure
        (or any earlier PlotParams) can leave them in the wrong palette. Call
        this immediately before building a figure to make the palette
        deterministic rather than dependent on import order.
        """
        if self.palatte == "dark":
            self.apply_dark()
        if self.palatte == "light":
            self.apply_light()

        self.apply_scaling()

        self.apply_grid()

        self.apply_font()
        return self


    def apply_font(self):
        """Applies this instance's font family, math included.

        Runs last, after the palette: `apply_light`/`apply_dark` reset the
        stylesheet, which would otherwise wipe these back to the matplotlib
        defaults.
        """
        plt.rcParams["font.family"] = self.font
        # Head the generic sans-serif list with it as well, so anything that
        # asks for the family rather than the name -- mathtext's upright
        # runs, other libraries' artists -- lands on the same face.
        plt.rcParams["font.sans-serif"] = [self.font, *self.FONT_FALLBACKS]
        plt.rcParams["mathtext.fontset"] = self.MATHTEXT_FONTSET
        return self


    def apply_grid(self):
        """Applies this instance's grid setting. Runs after the scaling so it
        is the single place the grid is decided, rather than each scaling
        method carrying its own copy."""
        plt.rcParams["axes.grid"] = self.grid
        if self.grid:
            plt.rcParams["grid.alpha"] = 0.5
            plt.rcParams["grid.linewidth"] = 0.5
        return self


    # Palette-dependent standards

    @property
    def file_suffix(self):
        """The palette tag appended to saved figure filenames (``_light`` /
        ``_dark``), so both palettes of the same figure coexist in a
        directory instead of overwriting each other."""
        return f"_{self.palatte}"

    @property
    def foreground(self):
        """The color for annotations drawn over a figure (scalebars, marker
        edges, overlay text): whichever of black/white contrasts with this
        palette's background."""
        return self._STANDARDS[self.palatte]["foreground"]

    @property
    def background(self):
        """This palette's figure/axes background color."""
        return self._STANDARDS[self.palatte]["background"]

    @property
    def bad_color(self):
        """The color a colormap should paint NaN/masked pixels. Matched to
        the background so unfitted spaxels read as empty sky rather than as
        a data value; without it matplotlib leaves them transparent, which
        only looks right until the figure is placed on another background."""
        return self._STANDARDS[self.palatte]["bad"]

    # Scaling-dependent standards
    #
    # rcParams cover the text matplotlib draws on its own (axis labels, tick
    # labels, titles). These expose the same table to text a call site draws
    # itself -- a colorbar label, a scalebar, a compass -- so those follow the
    # scaling too instead of carrying a hardcoded size.

    @property
    def label_size(self):
        """Point size for an axis label, colorbar label included."""
        return self._SCALING_FONTS[self.scaling]["axes_label"]

    @property
    def tick_size(self):
        """Point size for tick labels."""
        return self._SCALING_FONTS[self.scaling]["tick"]

    @property
    def annotation_size(self):
        """Point size for text drawn inside the axes: the scalebar label, the
        compass letters, anything captioning the data rather than the frame."""
        return self._SCALING_FONTS[self.scaling]["annotation"]

    def label_colorbar(self, colorbar, label):
        """Labels a colorbar at this instance's scaling.

        The label reads downward on a right-hand colorbar and is padded clear
        of the tick labels. The pad tracks the tick size because that is what
        it has to clear -- a pad fixed for 28pt ticks leaves a conspicuous gap
        beside 14pt ones.
        """
        colorbar.set_label(label, fontsize=self.label_size, rotation=270,
                           labelpad=1.8 * self.tick_size)
        return colorbar

    def colormap(self, name, lut=None):
        """Returns a named colormap with this palette's `bad_color` applied,
        the standard way to build a colormap for a spaxel map.

        Arguments
        ---------
        name : str
            Any matplotlib colormap name
        lut : int, optional
            Number of discrete levels to resample the colormap to, for maps
            of integer quantities (e.g. a component count)
        """
        cmap = plt.get_cmap(name, lut).copy() if lut else plt.get_cmap(name).copy()
        cmap.set_bad(self.bad_color)
        return cmap

    def colors(self):
        """The categorical color cycle for this palette."""
        return self.dark_colors() if self.palatte == "dark" else self.light_colors()


    def hide_axes_labels(self):
        plt.rcParams["axes.labelsize"] = 0
        plt.rcParams["xtick.labelsize"] = 0
        plt.rcParams["ytick.labelsize"] = 0
        plt.rcParams["axes.grid"] = False

    def apply_scaling(self):
        """Applies this instance's scaling to every text-size rcParam.

        Every size a figure draws is set here, `axes.titlesize` included. It
        used to be set only by the palette methods, which run *before* the
        scaling in `apply()` and so were never overridden -- a `paper` figure
        came out with 16pt axis labels under a 28pt title.
        """
        fonts = self._SCALING_FONTS[self.scaling]
        plt.rcParams["axes.labelsize"] = fonts["axes_label"]
        plt.rcParams["axes.titlesize"] = fonts["title"]

        plt.rcParams["xtick.labelsize"] = fonts["tick"]
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.labelsize"] = fonts["tick"]
        plt.rcParams["ytick.direction"] = "in"

        plt.rcParams["legend.fontsize"] = fonts["legend"]

        plt.rcParams["savefig.dpi"] = 1200
        plt.rcParams["savefig.bbox"] = "tight"
        return self


    def apply_light(self):
        """
        Applies light mode settings.

        Colors and canvas only: text sizes belong to the scaling, and
        `apply()` runs the scaling after the palette, so any size set here
        would silently win over the caller's choice of scaling.
        """

        # Reset to the stock stylesheet first. Without this, a light palette
        # applied after a dark one in the same process inherits every color
        # `dark_background` set (white text and axes on a black figure), and
        # the "light" figure comes out unreadable.
        plt.style.use("default")
        plt.rcParams["figure.figsize"] = (8,6)

        # Saved figures carry the palette's background rather than being
        # transparent, so a dark figure never renders as dark-on-dark (and a
        # light one never as light-on-light) in a viewer or slide deck.
        plt.rcParams["savefig.transparent"] = False
        plt.rcParams["savefig.facecolor"] = self.background
        plt.rcParams["savefig.edgecolor"] = self.background


    def apply_dark(self):
        """
        Applies dark mode settings. Colors and canvas only, as in
        `apply_light`.
        """

        # "default" first so the dark stylesheet lands on stock rcParams
        # rather than on whatever a previous palette left behind.
        plt.style.use(["default", "dark_background"])
        plt.rcParams["figure.figsize"] = (8,6)

        plt.rcParams["savefig.transparent"] = False
        plt.rcParams["savefig.facecolor"] = self.background
        plt.rcParams["savefig.edgecolor"] = self.background


    def dark_colors(self):
        return ["#e60049", "#0bb4ff", "#ffa300", "#50e991", "#9b19f5", "#e6d800", "#dc0ab4", "#b3d4ff", "#00bfa0"]

    def light_colors(self):
        return plt.get_cmap('tab20')

    def linestyles(self):
        return list(matplotlib.lines.lineStyles.keys())

    def markers(self):
        return list(matplotlib.markers.MarkerStyle.markers.keys())
