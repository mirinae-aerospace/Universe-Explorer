from tkinter import *
from tkinter import messagebox, filedialog
from tkinter import ttk
import sv_ttk
from PIL import Image, ImageTk
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

#메인 윈도우
window = Tk()
window.title("Universe Explorer")
sv_ttk.set_theme("dark")

titlef = ["Malgun Gothic", 30, "bold"]
subtitlef = ["Malgun Gothic", 17, "bold"]
textf = ["Malgun Gothic", 15]

style = ttk.Style()
style.configure('TButton', font=textf)

def toggle_theme():
    if sv_ttk.get_theme() == "dark":
        sv_ttk.use_light_theme()
    elif sv_ttk.get_theme() == "light":
        sv_ttk.use_dark_theme()

title = ttk.Label(window, text="Universe Explorer", font=titlef)
title.grid(row=0,column=1)
l1 = ttk.Label(window, text="Universe Explorer와 함께 우주를 탐험해보세요.\nv.alpha25.1.1 (KOR)", font=textf)
l1.grid(row=1,column=1)

#천체 검색
def object_search():
    os_window = Toplevel(window)
    os_window.title("Object Search")

    title = ttk.Label(os_window, text="천체 검색", font=titlef)
    title.grid(row=0,column=1)

    iq = ttk.Label(os_window, text="식별자 검색", font=subtitlef)
    iq.grid(row=1,column=1)

    l1 = ttk.Label(os_window, text="천체 이름", font=textf)
    l1.grid(row=2,column=0)

    e1 = ttk.Entry(os_window)
    e1.grid(row=2,column=1)

    cq = ttk.Label(os_window, text="좌표 검색", font=subtitlef)
    cq.grid(row=3,column=1)

    l2 = ttk.Label(os_window, text="ra", font=textf)
    l2.grid(row=4,column=0)

    e2 = ttk.Entry(os_window)
    e2.grid(row=4, column=1)

    l3 = ttk.Label(os_window, text="dec", font=textf)
    l3.grid(row=5,column=0)

    e3 = ttk.Entry(os_window)
    e3.grid(row=5, column=1)

    l5 = ttk.Label(os_window, text="반경 (°) (0 ~ 30)", font=textf)
    l5.grid(row=6,column=0)

    s1 = ttk.Scale(os_window, from_=0, to=30, length=180, orient=HORIZONTAL)
    s1.set(2)
    s1.grid(row=6,column=1)

    t1 = Text(os_window, width=50, height=17, font=textf)
    t1.grid(row=7, column=1)
    t1.insert("1.0", "식별자 검색과 좌표 검색은 독립적으로 작동합니다.\n여기에 검색 결과가 나타납니다.")

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

    #식별자 검색 Identifier Search
    def simbad_search():
        t1.delete("1.0", "end")
        t1.insert("1.0", "SIMBAD 검색 중...")
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
                        t1.insert("1.0", f"'{object_name}' 에 대한 SIMBAD 검색 결과가 없습니다")
                        os_window.update()
                        return
                    main_id = result_table['main_id'][0]
                    otype = result_table['otype'][0]
                    def get_otype_description(otype_code):
                        return simbad_otypes.get(otype_code, f"{otype_code} (설명 없음)")
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
                        f"천체 유형: {des}\n"
                        f"ra (ICRSd): {ra}\n"
                        f"dec (ICRSd): {dec}\n"
                        f"스펙트럼 유형: {sptype}\n"
                        f"은하 유형: {morphtype}\n"
                        f"겉보기 등급 (가시광선 V): --\n"
                        f"겉보기 등급 (B): --\n"
                        f"색지수 (B-V): --\n"
                        f"시선 속도 (km/s): {rv}\n"
                        f"적색편이: {z}\n"
                        f"시차 (mas): {plx}"
                        )
                    t1.insert("1.0", t1text)
                    SDSSImageApp(window, ra, dec)
                    return
                main_id = result_table['main_id'][0]
                otype = result_table['otype'][0]
                def get_otype_description(otype_code):
                    return simbad_otypes.get(otype_code, f"{otype_code} (설명 없음)")
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
                    f"천체 유형: {des}\n"
                    f"ra (ICRSd): {ra}\n"
                    f"dec (ICRSd): {dec}\n"
                    f"스펙트럼 유형: {sptype}\n"
                    f"은하 유형: {morphtype}\n"
                    f"겉보기 등급 (가시광선 V): --\n"
                    f"겉보기 등급 (B): {bmag}\n"
                    f"색지수 (B-V): --\n"
                    f"시선 속도 (km/s): {rv}\n"
                    f"적색편이: {z}\n"
                    f"시차 (mas): {plx}"
                    )
                t1.insert("1.0", t1text)
                SDSSImageApp(window, ra, dec)
                return
        main_id = result_table['main_id'][0]
        otype = result_table['otype'][0]
        def get_otype_description(otype_code):
            return simbad_otypes.get(otype_code, f"{otype_code} (설명 없음)")
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
            f"천체 유형: {des}\n"
            f"ra (ICRSd): {ra}\n"
            f"dec (ICRSd): {dec}\n"
            f"스펙트럼 유형: {sptype}\n"
            f"은하 유형: {morphtype}\n"
            f"겉보기 등급 (가시광선 V): {flux}\n"
            f"겉보기 등급 (B): {bmag}\n"
            f"색지수 (B-V): {bv_index}\n"
            f"시선 속도 (km/s): {rv}\n"
            f"적색편이: {z}\n"
            f"시차 (mas): {plx}"
            )
        t1.insert("1.0", t1text)
        SDSSImageApp(window, ra, dec)
    #좌표 검색 Coordinate Search
    def coords_search():
        t1.delete("1.0", "end")
        t1.insert("1.0", "SIMBAD 검색 중...")
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
        SDSSImageApp(window, ra, dec)

    b1 = ttk.Button(os_window, text="SIMBAD 검색", command=simbad_search)
    b1.grid(row=2,column=2)

    b2 = ttk.Button(os_window, text="좌표 검색", command=coords_search)
    b2.grid(row=6,column=2)

#SDSS 이미지 뷰어 SDSS Image Viewer
class SDSSImageApp:
    def __init__(self, master, ra, dec):
        self.root = Toplevel(master)
        self.root.title("SDSS Image Viewer")

        self.mainid_Label = ttk.Label(self.root, text=f"ra {ra}, \ndec {dec}", font=titlef)
        self.mainid_Label.grid(row=0,column=1)

        self.scale_label = ttk.Label(self.root, text="스케일 (0 ~ 15)", font=textf)
        self.scale_label.grid(row=0,column=0)

        self.scale_slider = ttk.Scale(self.root, from_=0, to=15, orient=HORIZONTAL)
        self.scale_slider.set(5.0)
        self.scale_slider.grid(row=1,column=0)

        self.fetch_btn = ttk.Button(self.root, text="이미지 불러오기", command=self.fetch_image)
        self.fetch_btn.grid(row=0,column=2)


        self.save_btn = ttk.Button(self.root, text="이미지 저장하기", command=self.save_image)
        self.save_btn.grid(row=1,column=2)
        # 이미지 출력
        self.image_label = ttk.Label(self.root, text="SDSS 이미지 불러오는 중...")
        self.image_label.grid(row=2,column=1)

        # 현재 이미지 저장용
        self.current_image = None
        self.tk_image = None

        self.root.after(100, self.fetch_image)

    def fetch_image(self):
        self.ra=ra
        self.dec=dec
        self.scale = float(self.scale_slider.get())
        try:
            # SDSS 컷아웃 API 요청
            url = (
                f"https://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg"
                f"?ra={self.ra}&dec={self.dec}&scale={self.scale}&width=512&height=512&opt=G"
            )
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                self.image_label.configure(text=f"HTTP Error: {response.status_code}")
                messagebox.showerror("Error", f"HTTP Error: {response.status_code}")
                return

            img_data = response.content
            img_pil = Image.open(io.BytesIO(img_data))
            self.current_image = img_pil.copy()

            img_tk = ImageTk.PhotoImage(img_pil)

            self.tk_image = img_tk

            self.image_label.configure(image=self.tk_image)
            self.image_label.image = self.tk_image

        except Exception as e:
            self.image_label.configure(text=f"이미지 로드 실패: {str(e)}")
            messagebox.showerror("Error", str(e))

    def save_image(self):
        if self.current_image is None:
            messagebox.showwarning("Warning", "이미지를 먼저 불러오세요.")
            return

        default_name = f"SDSS_{main_id}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.jpg"
        file_path = filedialog.asksaveasfilename(defaultextension=".jpg", initialfile=default_name,
                                                filetypes=[("JPEG files", "*.jpg"), ("All files", "*.*")])
        if file_path:
            try:
                self.current_image.save(file_path, format='JPEG')
                messagebox.showinfo("Saved", f"이미지가 저장되었습니다.:\n{file_path}")
            except Exception as e:
                messagebox.showerror("Save Error", str(e))

b1 = ttk.Button(window, text="천체 검색", command=object_search)
b1.grid(row=2,column=0)

#3차원 우주 지도
def spacemap():
    window = Toplevel()
    window.title("3D Space Map")

    title = ttk.Label(window, text="3D 우주 지도", font=titlef)
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

    b1 = ttk.Button(window,  text="SDSS 3D 우주 지도", command=sdssmap)
    b1.grid(row=6,column=1)

b2 = ttk.Button(window, text="3D 우주 지도", command=spacemap)
b2.grid(row=2,column=3)

#소개
def about():
    window = Toplevel()
    window.title("About")

    title = ttk.Label(window, text="프로그램 소개", font=titlef)
    title.grid(row=0,column=1)

    l1 = ttk.Label(window, text="Universe Explorer v.alpha25.1.1 (KOR)\nUniverse Explorer는 Python Tkinter로 제작된 우주를 탐험해 볼 수 있는 프로그램입니다. \n" \
    "SIMBAD와 SDSS 등 천문 데이터에 접근하고 시각적으로 보여줄 수 있습니다.\n" \
    "천체 검색 및 3D 우주 지도 등을 할 수 있습니다.\n" \
    "이 프로그램을 통해 태양계 밖 거의 모든 천체의 정보를 찾을 수 있고, 우주를 3차원으로 시각화할 수 있습니다.", font=textf)
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

    title = ttk.Label(window, text="버전 정보", font=titlef)
    title.grid(row=0,column=1)

    l1 = ttk.Label(window, text="Universe Explorer v.alpha25.1.1 (KOR)\n"
    "programming by jake (mirinae aerospace)", font=textf)
    l1.grid(row=1,column=1)

def whatsnew():
    window = Toplevel()
    window.title("What's new")

    title = ttk.Label(window, text="새로운 기능", font=titlef)
    title.grid(row=0,column=1)

    l1 = ttk.Label(window, text="Universe Explorer v.alpha25.1.1 (KOR)\n" \
    "테마 바꾸기 메뉴 추가", font=textf)
    l1.grid(row=1,column=1)

menu = Menu(window)
info = Menu(menu, tearoff=0)
info.add_command(label="소개", command=about)
info.add_command(label="버전 정보", command=vinfo)
info.add_command(label="새로운 기능", command=whatsnew)
info.add_command(label="테마 바꾸기", command=toggle_theme)
menu.add_cascade(label="정보...", menu=info)
window.config(menu=menu)

window.mainloop()
