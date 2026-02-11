from tkinter import *
from tkinter import messagebox, filedialog
import tkinter as tk
from tkinter import ttk
import sv_ttk
from PIL import Image, ImageTk, ImageEnhance
import io
import requests
import datetime
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.cosmology import Planck18 as cosmo
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astroquery.simbad import Simbad
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import lightkurve as lk

#메인 윈도우
window = Tk()
window.title("Universe Explorer")
sv_ttk.set_theme("dark")

titlef = ["Arial", 30, "bold"]
subtitlef = ["Arial", 17, "bold"]
textf = ["Arial", 15]

style = ttk.Style()
style.configure('TButton', font=textf)

title = ttk.Label(window, text="Universe Explorer", font=titlef)
title.grid(row=0,column=1)
l1 = ttk.Label(window, text="Explore The Universe with Universe Explorer.\nv.alpha26.1 (ENG)", font=textf)
l1.grid(row=1,column=1)

def toggle_theme():
    if sv_ttk.get_theme() == "dark":
        sv_ttk.use_light_theme()
    elif sv_ttk.get_theme() == "light":
        sv_ttk.use_dark_theme()

#천체 검색
def object_search():
    os_window = Toplevel(window)
    os_window.title("Object Search")
    os_window.geometry("757x477+0+0")

    title = ttk.Label(os_window, text="Object Search", font=titlef)
    title.grid(row=0,column=1)

    iq = ttk.Label(os_window, text="Identifier Search", font=subtitlef)
    iq.grid(row=1,column=1)

    l1 = ttk.Label(os_window, text="Object Name", font=textf)
    l1.grid(row=2,column=0)

    e1 = ttk.Entry(os_window)
    e1.grid(row=2,column=1)

    cq = ttk.Label(os_window, text="Coordinate Search", font=subtitlef)
    cq.grid(row=3,column=1)

    l2 = ttk.Label(os_window, text="RA", font=textf)
    l2.grid(row=4,column=0)

    e2 = ttk.Entry(os_window)
    e2.grid(row=4, column=1)

    l3 = ttk.Label(os_window, text="Dec", font=textf)
    l3.grid(row=5,column=0)

    e3 = ttk.Entry(os_window)
    e3.grid(row=5, column=1)

    l5 = ttk.Label(os_window, text="radius (°) (0 ~ 30)", font=textf)
    l5.grid(row=6,column=0)

    s1 = ttk.Scale(os_window, from_=0, to=30, length=180, orient=HORIZONTAL)
    s1.set(2)
    s1.grid(row=6,column=1)

    t1 = Text(os_window, width=50, height=17, font=textf)
    t1.grid(row=7, column=1)
    t1.insert("1.0", "Identifier Search and Coordinate Search operate independently.\nSearch results appear here.")

    global main_id
    global ra
    global dec

    simbad_otypes = {
        "*": "Star",
        "Ma*": "Massive Star",
        "bC*": "beta Cep Variable",
        "sg*": "Evolved Supergiant",
        "s*r": "Red Supergiant",
        "s*y": "Yellow Supergiant",
        "s*b": "Blue Supergiant",
        "WR*": "Wolf-Rayet",
        "N*": "Neutron Star",
        "Psr": "Pulsar",
        "Y*O": "Young Stellar Object",
        "TT*": "T Tauri Star",
        "Ae*": "Herbig Ae/Be Star",
        "out": "Outflow",
        "HH": "Herbig-Haro Object",
        "MS*": "Main Sequence Star",
        "Be*": "Be Star",
        "BS*": "Blue Straggler",
        "SX*": "SX Phe Variable",
        "gD*": "gamma Dor Variable",
        "dS*": "delta Sct Variable",
        "Ev*": "Evolved Star",
        "RG*": "Red Giant Branch star",
        "HS*": "Hot Subdwarf",
        "HB*": "Horizontal Branch Star",
        "RR*": "RR Lyrae Variable",
        "WV*": "Type II Cepheid Variable",
        "Ce*": "Cepheid Variable",
        "cC*": "Classical Cepheid Variable",
        "C*": "Carbon Star",
        "S*": "S Star",
        "LP*": "Long-Period Variable",
        "AB*": "Asymptotic Giant Branch Star",
        "Mi*": "Mira Variable",
        "OH*": "OH/IR Star",
        "pA*": "Post-AGB Star",
        "RV*": "RV Tauri Variable",
        "PN": "Planetary Nebula",
        "WD*": "White Dwarf",
        "Pe*": "Chemically Peculiar Star",
        "a2*": "alpha2 CVn Variable",
        "RC*": "R CrB Variable",
        "**": "Double or Multiple Star",
        "EB*": "Eclipsing Binary",
        "El*": "Ellipsoidal Variable",
        "SB*": "Spectroscopic Binary",
        "RS*": "RS CVn Variable",
        "BY*": "BY Dra Variable",
        "Sy*": "Symbiotic Star",
        "XB*": "X-ray Binary",
        "LXB": "Low Mass X-ray Binary",
        "HXB": "High Mass X-ray Binary",
        "CV*": "Cataclysmic Binary",
        "No*": "Classical Nova",
        "SN*": "SuperNova",
        "LM*": "Low-mass Star",
        "BD*": "Brown Dwarf",
        "Pl": "Extra-solar Planet",
        "V*": "Variable Star",
        "Ir*": "Irregular Variable",
        "Er*": "Eruptive Variable",
        "Ro*": "Rotating Variable",
        "Pu*": "Pulsating Variable",
        "Em*": "Emission-line Star",
        "PM*": "High Proper Motion Star",
        "HV*": "High Velocity Star",
        "Cl*": "Cluster of Stars",
        "GlC": "Globular Cluster",
        "OpC": "Open Cluster",
        "As*": "Association of Stars",
        "St*": "Stellar Stream",
        "MGr": "Moving Group",
        "ISM": "Interstellar Medium Object",
        "SFR": "Star Forming Region",
        "HII": "HII Region",
        "Cld": "Cloud",
        "GNe": "Nebula",
        "RNe": "Reflection Nebula",
        "MoC": "Molecular Cloud",
        "DNe": "Dark Cloud",
        "glb": "Globule",
        "CGb": "Cometary Globule",
        "HVC": "High-velocity Cloud",
        "cor": "Dense Core",
        "bub": "Bubble",
        "SNR": "SuperNova Remnant",
        "sh": "Interstellar Shell",
        "flt": "Interstellar Filament",
        "G": "Galaxy",
        "LSB": "Low Surface Brightness Galaxy",
        "bCG": "Blue Compact Galaxy",
        "SBG": "Starburst Galaxy",
        "H2G": "HII Galaxy",
        "EmG": "Emission-line galaxy",
        "AGN": "Active Galaxy Nucleus",
        "SyG": "Seyfert Galaxy",
        "Sy1": "Seyfert 1 Galaxy",
        "Sy2": "Seyfert 2 Galaxy",
        "rG": "Radio Galaxy",
        "LIN": "LINER-type Active Galaxy Nucleus",
        "QSO": "Quasar",
        "Bla": "Blazar",
        "BLL": "BL Lac",
        "GiP": "Galaxy in Pair of Galaxies",
        "GiG": "Galaxy towards a Group of Galaxies",
        "GiC": "Galaxy towards a Cluster of Galaxies",
        "BiC": "Brightest Galaxy in a Cluster",
        "IG": "Interacting Galaxies",
        "PaG": "Pair of Galaxies",
        "GrG": "Group of Galaxies",
        "CGG": "Compact Group of Galaxies",
        "ClG": "Cluster of Galaxies",
        "PCG": "Proto Cluster of Galaxies",
        "SCG": "Supercluster of Galaxies",
        "vid": "Underdense Region of the Universe",
        "grv": "Gravitational Source",
        "Lev": "(Micro)Lensing Event",
        "gLS": "Gravitational Lens System",
        "gLe": "Gravitational Lens",
        "LeI": "Gravitationally Lensed Image",
        "LeG": "Gravitationally Lensed Image of a Galaxy",
        "LeQ": "Gravitationally Lensed Image of a Quasar",
        "BH": "Black Hole",
        "GWE": "Gravitational Wave Event",
        "ev": "Transient Event",
        "var": "Variable source",
        "Rad": "Radio Source",
        "mR": "Metric Radio Source",
        "cm": "Centimetric Radio Source",
        "mm": "Millimetric Radio Source",
        "smm": "Sub-Millimetric Source",
        "HI": "HI (21cm) Source",
        "rB": "Radio Burst",
        "Mas": "Maser",
        "IR": "Infra-Red Source",
        "FIR": "Far-IR source",
        "MIR": "Mid-IR Source",
        "NIR": "Near-IR Source",
        "Opt": "Optical Source",
        "EmO": "Emission Object",
        "blu": "Blue Object",
        "UV": "UV-emission Source",
        "X": "X-ray Source",
        "ULX": "Ultra-luminous X-ray Source",
        "gam": "Gamma-ray Source",
        "gB": "Gamma-ray Burst",
        "mul": "Composite Object, Blend",
        "err": "Not an Object (Error, Artefact, ...)",
        "PoC": "Part of Cloud",
        "PoG": "Part of a Galaxy",
        "?": "Object of Unknown Nature",
        "reg": "Region defined in the Sky"
    }

    def results():
        TwoMASSImageApp(os_window, ra, dec)
        re_win = Toplevel()
        re_win.title("Search Results")
        l1 = ttk.Label(re_win, text="Search Results", font=titlef)
        l1.grid(row=0,column=1)
        b1 = ttk.Button(re_win, text="View 2MASS Image", command=lambda: TwoMASSImageApp(re_win, ra, dec))
        b1.grid(row=1,column=0)
        b2 = ttk.Button(re_win, text="View Kepler Light Curve", command=lambda: LightCurveApp(Toplevel(), main_id))
        b2.grid(row=1,column=1)
        b3 = ttk.Button(re_win, text="View Magnitudes(planned)", command=lambda: MagnitudeApp(re_win, main_id))
        b3.grid(row=1,column=2)

    #식별자 검색 Identifier Search
    def simbad_search():
        t1.delete("1.0", "end")
        t1.insert("1.0", "Searching SIMBAD...")
        os_window.update()
        object_name = str(e1.get())
        simbad = Simbad()
        result_table = None
        t1text = None
        simbad.reset_votable_fields()
        simbad.add_votable_fields('otype', 'sp_type', 'morph_type', 'V', 'B', 'rvz_radvel', 'rvz_redshift', 'plx_value')
        result_table = simbad.query_object(object_name)
        global main_id
        global ra
        global dec
        #검색 결과 없을 때 When No Results
        if result_table is None or len(result_table) == 0:
                t1.delete("1.0", "end")
                t1.insert("1.0", f"NoResultsWarning: The request executed correctly, \nbut there was no data corresponding to these criteria in SIMBAD: '{object_name}'. Retrying...")
                os_window.update()
                simbad.reset_votable_fields()
                simbad.add_votable_fields('otype', 'sp_type', 'B', 'morph_type', 'rvz_radvel', 'rvz_redshift', 'plx_value')
                result_table = simbad.query_object(object_name)
                if result_table is None or len(result_table) == 0:
                    t1.delete("1.0", "end")
                    t1.insert("1.0", f"NoResultsWarning: The request executed correctly, \nbut there was no data corresponding to these criteria in SIMBAD: '{object_name}'. Retrying2...")
                    os_window.update()
                    simbad.reset_votable_fields()
                    simbad.add_votable_fields('otype', 'sp_type', 'morph_type', 'rvz_radvel', 'rvz_redshift', 'plx_value')
                    result_table = simbad.query_object(object_name)
                    if result_table is None or len(result_table) == 0:
                        t1.delete("1.0", "end")
                        t1.insert("1.0", f"There is No Result in SIMBAD '{object_name}'")
                        os_window.update()
                        return
                    main_id = result_table['main_id'][0]
                    otype = result_table['otype'][0]
                    def get_otype_description(otype_code):
                        return simbad_otypes.get(otype_code, f"{otype_code} (No Description)")
                    des = get_otype_description(otype)
                    ra = result_table['ra'][0]
                    dec = result_table['dec'][0]
                    sptype = result_table['sp_type'][0]
                    morphtype = result_table['morph_type'][0]
                    rv = result_table['rvz_radvel'][0]
                    z = result_table['rvz_redshift'][0]
                    plx = result_table['plx_value'][0]
                    t1.delete("1.0", "end")
                    t1text = (
                        f"{main_id}\n"
                        f"Object Type: {des}\n"
                        f"ra (ICRSd): {ra}\n"
                        f"dec (ICRSd): {dec}\n"
                        f"Spectrum Type: {sptype}\n"
                        f"Morphological Type: {morphtype}\n"
                        f"Apparent Magnitude (V): --\n"
                        f"Apparent Magnitude (B): --\n"
                        f"Color index (B-V): --\n"
                        f"Radial velocity (km/s): {rv}\n"
                        f"Redshift: {z}\n"
                        f"Parallax (mas): {plx}"
                        )
                    t1.insert("1.0", t1text)
                    results()
                    return
                main_id = result_table['main_id'][0]
                otype = result_table['otype'][0]
                def get_otype_description(otype_code):
                    return simbad_otypes.get(otype_code, f"{otype_code} (No Description)")
                des = get_otype_description(otype)
                ra = result_table['ra'][0]
                dec = result_table['dec'][0]
                sptype = result_table['sp_type'][0]
                morphtype = result_table['morph_type'][0]
                bmag = result_table['B'][0]
                rv = result_table['rvz_radvel'][0]
                z = result_table['rvz_redshift'][0]
                plx = result_table['plx_value'][0]
                t1.delete("1.0", "end")
                t1text = (
                    f"{main_id}\n"
                    f"Object Type: {des}\n"
                    f"ra (ICRSd): {ra}\n"
                    f"dec (ICRSd): {dec}\n"
                    f"Spectrum Type: {sptype}\n"
                    f"Morphological Type: {morphtype}\n"
                    f"Apparent Magnitude (V): --\n"
                    f"Apparent Magnitude (B): {bmag}\n"
                    f"Color index (B-V): --\n"
                    f"Radial velocity (km/s): {rv}\n"
                    f"Redshift: {z}\n"
                    f"Parallax (mas): {plx}"
                    )
                t1.insert("1.0", t1text)
                results()     
                return
        main_id = result_table['main_id'][0]
        otype = result_table['otype'][0]
        def get_otype_description(otype_code):
            return simbad_otypes.get(otype_code, f"{otype_code} (No Description)")
        des = get_otype_description(otype)
        ra = result_table['ra'][0]
        dec = result_table['dec'][0]
        sptype = result_table['sp_type'][0]
        morphtype = result_table['morph_type'][0]
        flux = result_table['V'][0]
        bmag = result_table['B'][0]
        bv_index = "N/A"
        try:
            if flux != "N/A" and bmag != "N/A":
                bv_value = float(bmag) - float(flux)
                bv_index = f"{bv_value:.3f}"
        except ValueError:
            pass
        rv = result_table['rvz_radvel'][0]
        z = result_table['rvz_redshift'][0]
        plx = result_table['plx_value'][0]
        t1.delete("1.0", "end")
        t1text = (
            f"{main_id}\n"
            f"Object Type: {des}\n"
            f"ra (ICRSd): {ra}\n"
            f"dec (ICRSd): {dec}\n"
            f"Spectrum Type: {sptype}\n"
            f"Morphological Type: {morphtype}\n"
            f"Apparent Magnitude (V): {flux}\n"
            f"Apparent Magnitude (B): {bmag}\n"
            f"Color index (B-V): {bv_index}\n"
            f"Radial velocity (km/s): {rv}\n"
            f"Redshift: {z}\n"
            f"Parallax (mas): {plx}"
            )
        t1.insert("1.0", t1text)
        results()
    #좌표 검색 Coordinate Search
    def coords_search():
        t1.delete("1.0", "end")
        t1.insert("1.0", "Searching SIMBAD...")
        os_window.update()
        global ra
        global dec
        ra = float(e2.get())
        dec = float(e3.get())
        radius = float(s1.get())
        simbad = Simbad()
        simbad.ROW_LIMIT = 30000
        simbad.add_votable_fields("otype")
        result_table = simbad.query_region(SkyCoord(ra, dec, unit=(u.deg, u.deg)), radius=radius * u.arcsec) #, criteria="otype = 'Galaxy..'"
        t1text = result_table[["main_id", "ra", "dec"]]
        t1.delete("1.0", "end")
        t1.insert("1.0", t1text)
        TwoMASSImageApp(os_window, ra, dec)

    b1 = ttk.Button(os_window, text="Search SIMBAD", command=simbad_search)
    b1.grid(row=2,column=2)

    b2 = ttk.Button(os_window, text="Coordinates Search", command=coords_search)
    b2.grid(row=6,column=2)

#2MASS 흑백 이미지 뷰어
class TwoMASSImageApp:
    def __init__(self, master, ra, dec):
        self.iv_win = Toplevel(master)
        self.iv_win.title(f"2MASS Image Viewer: {main_id}")
        self.iv_win.geometry("490x630+757+0")
        
        self.iv_win.columnconfigure(1, weight=1)
        self.iv_win.rowconfigure(4, weight=1)

        self.ra, self.dec = ra, dec
        self.raw_pil_img = None 

        ttk.Label(self.iv_win, text=f"RA: {ra},\nDec: {dec}", font=titlef).grid(row=0, column=0, columnspan=4, pady=5)

        ttk.Label(self.iv_win, text="Image range:").grid(row=1, column=0, padx=5, sticky=E)
        self.size_slider = ttk.Scale(self.iv_win, from_=1, to=40, orient=HORIZONTAL)
        self.size_slider.set(15.0)
        self.size_slider.grid(row=1, column=1, padx=5, sticky=EW)

        ttk.Label(self.iv_win, text="Brightness:").grid(row=2, column=0, padx=5, sticky=E)
        self.bright_slider = ttk.Scale(self.iv_win, from_=0.1, to=3.0, orient=HORIZONTAL, command=lambda e: self.apply_filter())
        self.bright_slider.set(1.0)
        self.bright_slider.grid(row=2, column=1, padx=5, sticky=EW)

        ttk.Label(self.iv_win, text="Contrast:").grid(row=3, column=0, padx=5, sticky=E)
        self.contrast_slider = ttk.Scale(self.iv_win, from_=0.5, to=5.0, orient=HORIZONTAL, command=lambda e: self.apply_filter())
        self.contrast_slider.set(1.0)
        self.contrast_slider.grid(row=3, column=1, padx=5, sticky=EW)

        btn_frame = ttk.Frame(self.iv_win)
        btn_frame.grid(row=1, column=2, rowspan=3, padx=10)
        ttk.Button(btn_frame, text="Get Image", command=self.fetch_image).pack(fill=X)
        ttk.Button(btn_frame, text="Save Image", command=self.save_image).pack(fill=X, pady=5)

        self.img_label = ttk.Label(self.iv_win, text="Waiting...", background="black", foreground="white")
        self.img_label.grid(row=4, column=0, columnspan=4, padx=10, pady=10, sticky=NSEW)

        self.iv_win.after(100, self.fetch_image)

    def fetch_image(self):
        try:
            pos_query = f"{self.ra} {self.dec}".strip()
            deg = self.size_slider.get() / 60.0
            url = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
            params = {"survey": "2massk", "position": pos_query, "size": deg, "pixels": 600, "return": "jpg"}

            resp = requests.get(url, params=params, timeout=15)
            if resp.status_code == 200:
                self.raw_pil_img = Image.open(io.BytesIO(resp.content)).convert("L")
                self.apply_filter()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def apply_filter(self):
        if self.raw_pil_img is None: return
        
        img = self.raw_pil_img.copy()
        
        enhancer_b = ImageEnhance.Brightness(img)
        img = enhancer_b.enhance(self.bright_slider.get())
        
        enhancer_c = ImageEnhance.Contrast(img)
        img = enhancer_c.enhance(self.contrast_slider.get())
        
        self.current_image = img
        
        preview = img.resize((480, 480), Image.Resampling.LANCZOS)
        self.tk_img = ImageTk.PhotoImage(preview)
        self.img_label.config(image=self.tk_img, text="")

    def save_image(self):
        if not self.current_image: return
        path = filedialog.asksaveasfilename(defaultextension=".jpg", initialfile=f"2MASS_{main_id}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.jpg")
        if path: self.current_image.save(path)

#케플러 광도곡선 뷰어
class LightCurveApp:
    def __init__(self, lc_win, star_id):
        self.lc_win = lc_win
        self.star_id = star_id
        self.lc_win.title(f"Kepler Light Curve Viewer '{self.star_id}'")

        screen_width = self.lc_win.winfo_screenwidth()
        screen_height = self.lc_win.winfo_screenheight()

        width = screen_width // 2
        height = int(screen_height * 0.9)

        x_pos = screen_width // 2
        y_pos = 0

        self.lc_win.geometry(f"{width}x{height}+{x_pos}+{y_pos}")

        self.style = ttk.Style()

        self.status_frame = ttk.Frame(self.lc_win)
        self.status_frame.pack(side=tk.TOP, fill=tk.X, padx=15, pady=5)

        self.status_label = ttk.Label(
            self.status_frame, 
            text=f"Waiting for analysis: {self.star_id}", 
            foreground="blue"
        )
        self.status_label.pack(side=tk.LEFT)

        if sv_ttk.get_theme() == "dark":
            self.fig, self.axes = plt.subplots(3, 1, figsize=(8, 10), facecolor="#303030")
        elif sv_ttk.get_theme() == "light":
            self.fig, self.axes = plt.subplots(3, 1, figsize=(8, 10), facecolor="#FFFFFF")
        plt.tight_layout(pad=5.0)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.lc_win)

        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.lc_win.after(100, self.process_data)

    def process_data(self):
        self.status_label.config(text=f"Loading and analyzing data... ({self.star_id})", foreground="red")
        self.lc_win.update_idletasks() 

        try:
            search_result = lk.search_targetpixelfile(self.star_id, quarter=4, author="Kepler")
            
            if len(search_result) == 0:
                self.status_label.config(text="No data found for target", foreground="black")
                messagebox.showwarning("Warning", f"No data found for target {self.star_id}")
                return

            tpf = search_result.download()
            lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)

            flat_lc = lc.flatten(window_length=401)
            periods = np.linspace(1.0, 10.0, 5000)
            periodogram = flat_lc.to_periodogram(method='bls', period=periods)
            best_period = periodogram.period_at_max_power
            folded_lc = flat_lc.fold(period=best_period)

            titles = [f"{self.star_id} Original", "Flattened", f"Folded (Period: {best_period:.4f} d)"]
            data_list = [lc, flat_lc, folded_lc]

            for ax, data, title in zip(self.axes, data_list, titles):
                ax.clear()
                ax.set_facecolor("#FFFFFF")
                data.scatter(ax=ax, title=title)
                y_min = np.nanpercentile(data.flux.value, 0.5)
                y_max = np.nanpercentile(data.flux.value, 99.5)
                ax.set_ylim(y_min, y_max)

            self.canvas.draw()
            self.status_label.config(text=f"Analysis success: {self.star_id}", foreground="green")
            
        except Exception as e:
            self.status_label.config(text="Analysis error occurred", foreground="black")
            messagebox.showerror("Error", f"Analysis error: {e}")

class MagnitudeApp:
    pass

b1 = ttk.Button(window, text="Object Search", command=object_search)
b1.grid(row=2,column=0)

#3차원 우주 지도
def spacemap():
    window = Toplevel()
    window.title("3D Space Map")

    title = ttk.Label(window, text="3D Space Map", font=titlef)
    title.grid(row=0,column=1)

    l1 = ttk.Label(window, text="ra range between", font=textf)
    l1.grid(row=1,column=0)
    e1 = ttk.Entry(window)
    e1.grid(row=1,column=1)
    l2 = ttk.Label(window, text="and", font=textf)
    l2.grid(row=1,column=2)
    e2 = ttk.Entry(window)
    e2.grid(row=1,column=3)

    l3 = ttk.Label(window, text="dec range between", font=textf)
    l3.grid(row=2,column=0)
    e3 = ttk.Entry(window)
    e3.grid(row=2,column=1)
    l4 = ttk.Label(window, text="and", font=textf)
    l4.grid(row=2,column=2)
    e4 = ttk.Entry(window)
    e4.grid(row=2,column=3)

    l5 = ttk.Label(window, text="max results (0 ~ 10000)")
    l5.grid(row=3,column=0)
    s1 = ttk.Scale(window, from_=0, to=10000, length=180, orient=HORIZONTAL)
    s1.set(5000)
    s1.grid(row=3,column=1)

    l6 = ttk.Label(window, text="galaxy max z (0 ~ 15)")
    l6.grid(row=4,column=0)
    s2 = ttk.Scale(window, from_=0, to=15, length=180, orient=HORIZONTAL)
    s2.set(1)
    s2.grid(row=4,column=1)

    l7 = ttk.Label(window, text="qso max z")
    l7.grid(row=5,column=0)
    e5 = ttk.Entry(window)
    e5.grid(row=5,column=1)

    e1.insert(1,'0')
    e2.insert(1,'360')
    e3.insert(1,'-10')
    e4.insert(1,'70')
    e5.insert(1,'0.01')

    def sdssmap():
        def fetch_objects(object_class, ra_range, dec_range, max_results, min_z, max_z):
            query = f"""
            SELECT TOP {max_results} ra, dec, class, z
            FROM SpecObj
            WHERE class = '{object_class}'
            AND z BETWEEN {min_z} AND {max_z}
            AND ra BETWEEN {ra_range[0]} AND {ra_range[1]}
            AND dec BETWEEN {dec_range[0]} AND {dec_range[1]}
            """
            return SDSS.query_sql(query)

        def convert_to_xyz(ra, dec, z):
            dist_ly = cosmo.luminosity_distance(z).to(u.lyr).value
            ra_rad = np.deg2rad(ra)
            dec_rad = np.deg2rad(dec)
            x = dist_ly * np.cos(dec_rad) * np.cos(ra_rad)
            y = dist_ly * np.cos(dec_rad) * np.sin(ra_rad)
            z_coord = dist_ly * np.sin(dec_rad)
            return x, y, z_coord

        def plot_3d_map(galaxies, qsos):
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111, projection='3d')

            if len(galaxies) > 0:
                x_gal, y_gal, z_gal = convert_to_xyz(galaxies['ra'], galaxies['dec'], galaxies['z'])
                ax.scatter(x_gal, y_gal, z_gal, c="#0000FF", s=6, alpha=0.7, label='Galaxy', linewidths=0.3)

            if len(qsos) > 0:
                x_qso, y_qso, z_qso = convert_to_xyz(qsos['ra'], qsos['dec'], qsos['z'])
                ax.scatter(x_qso, y_qso, z_qso, c="#FF0000", s=5, alpha=0.6, label='QSO', linewidths=0.2)

            ax.set_xlabel('X [ly]', color='white')
            ax.set_ylabel('Y [ly]', color='white')
            ax.set_zlabel('Z [ly]', color='white')

            ax.set_title('SDSS 3D Space Map \nblue=galaxy, red=qso\nleft click for turn, right click for zoom in/out, scroll button for move', color='white')

            ax.tick_params(colors='white')

            ax.grid(False)

            fig.patch.set_facecolor('black')
            ax.set_facecolor('black')

            ax.set_box_aspect([1, 1, 1])

            plt.show()

        ra_range0 = float(e1.get())
        ra_range1 = float(e2.get())
        dec_range0 = float(e3.get())
        dec_range1 = float(e4.get())
        max_results = int(s1.get())
        galaxy_max_z = int(s2.get())
        qso_max_z = float(e5.get())

        galaxies = fetch_objects('GALAXY', ra_range=(ra_range0, ra_range1), dec_range=(dec_range0, dec_range1), max_results=max_results, min_z=0, max_z=galaxy_max_z)
        qsos = fetch_objects('QSO', ra_range=(ra_range0, ra_range1), dec_range=(dec_range0, dec_range1), max_results=max_results, min_z=0, max_z=qso_max_z)

        plot_3d_map(galaxies, qsos)

    b1 = ttk.Button(window,  text="SDSS 3D Space Map", command=sdssmap)
    b1.grid(row=6,column=1)

b2 = ttk.Button(window, text="3D Space Map", command=spacemap)
b2.grid(row=2,column=3)

#소개
def about():
    window = Toplevel()
    window.title("About")

    title = ttk.Label(window, text="About This Program", font=titlef)
    title.grid(row=0,column=1)

    l1 = ttk.Label(window, text="Universe Explorer v.alpha26.1 (ENG)\nUniverse Explorer is a program for exploring the universe, written in Python Tkinter.\n" \
    "It provides access to and visualization of astronomical data, including SIMBAD and SDSS.\n" \
    "It offers celestial object searches, and 3D space map.\n" \
    "You can also check photos of celestial objects and their light curves if they were observed by the Kepler Space Telescope.\n" \
    "This program allows you to find information on virtually every celestial object outside our solar system, \n" \
    "and visualize the universe in 3D.", font=textf)
    l1.grid(row=1,column=1)

    try:
        img = Image.open("about.gif")
        photo = ImageTk.PhotoImage(img)
        plabel = ttk.Label(window, image=photo)
        plabel.image = photo
        plabel.grid(row=2,column=1)
    except FileNotFoundError:
        plabel = ttk.Label(window, text="[Image Not Found: about.gif]", font=textf)
        plabel.grid(row=2,column=1)


#버전 정보
def vinfo():
    window = Toplevel()
    window.title("Version Infomation")

    title = ttk.Label(window, text="Version Infomation", font=titlef)
    title.grid(row=0,column=1)

    l1 = ttk.Label(window, text="Universe Explorer v.alpha26.1 (ENG)\n"
    "programming by jake (mirinae aerospace)", font=textf)
    l1.grid(row=1,column=1)

def whatsnew():
    window = Toplevel()
    window.title("What's new")

    title = ttk.Label(window, text="What's new", font=titlef)
    title.grid(row=0,column=1)

    l1 = ttk.Label(window, text="Universe Explorer v.alpha26.1 (ENG)\n" \
    "Added Kepler Light Curve Viewer, changed Survey in Image Viewer to 2MASS, and changed font.", font=textf)
    l1.grid(row=1,column=1)

menu = Menu(window)
info = Menu(menu, tearoff=0)
info.add_command(label="About", command=about)
info.add_command(label="Version Info", command=vinfo)
info.add_command(label="What's new", command=whatsnew)
info.add_command(label="Toggle Theme", command=toggle_theme)
menu.add_cascade(label="Info...", menu=info)
window.config(menu=menu)

window.mainloop()
