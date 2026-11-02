"""
Enriches `test files/chesapeake_bay_bridge_tunnel_4x3.ifc` with:
1. Wheeled Vehicles (`IfcVehicle`, VEHICLEWHEELED):
   - Passenger Sedans & Crossover SUVs (4 distinct multi-part types/colors) along NB & SB spans, approaches, toll plazas, portal islands, open-cut ramps, and tunnels
   - 18-Wheeler Tractor-Trailer Freight Trucks (`US-13 Delmarva Reefer & Container Semi-Trailers`)
   - CBBT Authority Emergency Patrol & Tunnel Response Tow Trucks parked on the paved pull-off shoulder aprons at the 4 Portal Islands and 2 Toll Plazas
2. Marine Vessels (`IfcVehicle`, VEHICLEMARINE):
   - Panamax Container Ship (`MV Chesapeake Bay Express`, 210m x 32m) sailing through the Thimble Shoal Shipping Channel above the Thimble Shoal Tunnel
   - Naval Auxiliary Cargo Vessel (`USNS Hampton Roads`, 175m x 28m) transiting the Chesapeake Navigation Channel above the Chesapeake Channel Tunnel
   - Channel Escort Tugboat (`MV Cape Charles`, 32m x 9.5m) passing beneath the North Channel High-Level Truss Bridge
   - Coastal Pilot & Charter Vessel (`MV Sea Gull II`, 22m x 6.8m) berthed alongside the Sea Gull Fishing Pier at Island #1
3. ITS Lane-Control Signals, Overheight Vehicle Laser Detectors & PTZ CCTV Cameras (`IfcSignal`, `IfcSensor`, `IfcAudioVisualAppliance`)
4. Tunnel Emergency & Life-Safety Systems (`IfcFireSuppressionTerminal`, `IfcAlarm`, `IfcUnitaryEquipment`)
"""

import math
import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.util.placement as pl
import numpy as np

IFC_PATH = "test files/chesapeake_bay_bridge_tunnel_4x3.ifc"


def new_guid() -> str:
    return ifcopenshell.guid.new()


def main():
    f = ifcopenshell.open(IFC_PATH)

    # Guard against duplicate runs
    existing_vehicles = f.by_type("IfcVehicle")
    if existing_vehicles:
        print(f"File already contains {len(existing_vehicles)} IfcVehicle entities; skipping duplicate insertion.")
        return

    owner_history = f.by_type("IfcOwnerHistory")[0]
    body_ctx = [c for c in f.by_type("IfcGeometricRepresentationSubContext") if c.ContextIdentifier == "Body"][0]
    origin_2d = f.createIfcAxis2Placement2D(
        f.createIfcCartesianPoint((0.0, 0.0)),
        f.createIfcDirection((1.0, 0.0)),
    )
    z_dir = f.createIfcDirection((0.0, 0.0, 1.0))
    origin_3d = f.createIfcAxis2Placement3D(
        f.createIfcCartesianPoint((0.0, 0.0, 0.0)),
        z_dir,
        f.createIfcDirection((1.0, 0.0, 0.0)),
    )
    ident_xform = f.by_type("IfcCartesianTransformationOperator3D")[0]

    # Helper to create or reuse surface style by name
    existing_styles = {s.Name: s for s in f.by_type("IfcSurfaceStyle") if getattr(s, "Name", None)}

    def make_style(name: str, r: float, g: float, b: float):
        if name in existing_styles:
            return existing_styles[name]
        col = f.createIfcColourRgb(name, float(r), float(g), float(b))
        sh = f.createIfcSurfaceStyleRendering(
            col,
            0.0,
            f.createIfcNormalisedRatioMeasure(0.65),
            None,
            None,
            None,
            None,
            None,
            "NOTDEFINED",
        )
        st = f.createIfcSurfaceStyle(name, "BOTH", [sh])
        existing_styles[name] = st
        return st

    style_navy_hull = make_style("Container Ship Deep Navy Hull", 0.10, 0.20, 0.36)
    style_antifouling = make_style("Marine Antifouling Boot-Top Red", 0.68, 0.16, 0.14)
    style_deck_white = make_style("Marine Superstructure Bridge White", 0.94, 0.95, 0.97)
    style_container_orange = make_style("Intermodal Container Hapag Orange", 0.90, 0.42, 0.10)
    style_container_teal = make_style("Intermodal Container Maersk Blue", 0.20, 0.58, 0.76)
    style_naval_grey = make_style("US Navy Haze Grey Hull & Superstructure", 0.56, 0.60, 0.64)
    style_car_blue = make_style("Automotive Atlantic Metallic Blue", 0.14, 0.38, 0.74)
    style_car_red = make_style("Automotive Crimson Red", 0.80, 0.14, 0.16)
    style_car_silver = make_style("Automotive Titanium Silver", 0.76, 0.79, 0.82)
    style_car_white = make_style("Automotive Pearl White", 0.93, 0.94, 0.96)
    style_cabin_glass = make_style("Automotive Tinted Windshield Glass", 0.18, 0.26, 0.34)
    style_tire_black = make_style("High-Energy Marine Elastomer Black", 0.15, 0.16, 0.18)
    style_safety_yellow = make_style("MUTCD Highway Centerline Yellow", 0.95, 0.78, 0.12)
    style_safety_orange = make_style("MUTCD Traffic Cone Safety Orange", 0.96, 0.40, 0.08)
    style_signal_red = make_style("ITS Lane Control Signal Red X", 0.90, 0.14, 0.14)
    style_signal_green = make_style("ITS Lane Control Signal Green Arrow", 0.12, 0.82, 0.32)
    style_galv_steel = make_style("Hot-Dip Galvanized Highway Steel", 0.70, 0.73, 0.76)
    style_sos_blue = make_style("CBBT Tunnel Portal Fascia Sign Blue", 0.10, 0.34, 0.68)
    style_jet_blue = make_style("Tunnel Jet Fan Ventilation Blue", 0.18, 0.45, 0.76)

    def make_box_solid(dx: float, dy: float, dz: float, cx: float, cy: float, cz: float, style, label: str):
        prof = f.createIfcRectangleProfileDef("AREA", label, origin_2d, float(dx), float(dy))
        pos = f.createIfcAxis2Placement3D(
            f.createIfcCartesianPoint((float(cx), float(cy), float(cz))),
            z_dir,
            f.createIfcDirection((1.0, 0.0, 0.0)),
        )
        solid = f.createIfcExtrudedAreaSolid(prof, pos, z_dir, float(dz))
        f.createIfcStyledItem(solid, [style], label)
        return solid

    def make_poly_solid(pts_2d, dz: float, cz: float, style, label: str):
        ifc_pts = [f.createIfcCartesianPoint((float(px), float(py))) for px, py in pts_2d]
        ifc_pts.append(ifc_pts[0])
        poly = f.createIfcPolyline(ifc_pts)
        prof = f.createIfcArbitraryClosedProfileDef("AREA", label, poly)
        pos = f.createIfcAxis2Placement3D(
            f.createIfcCartesianPoint((0.0, 0.0, float(cz))),
            z_dir,
            f.createIfcDirection((1.0, 0.0, 0.0)),
        )
        solid = f.createIfcExtrudedAreaSolid(prof, pos, z_dir, float(dz))
        f.createIfcStyledItem(solid, [style], label)
        return solid

    def make_rep_map(items):
        shape_rep = f.createIfcShapeRepresentation(body_ctx, "Body", "SweptSolid", items)
        return f.createIfcRepresentationMap(origin_3d, shape_rep)

    # -------------------------------------------------------------------------
    # 1. Build Multi-Part RepresentationMaps & ElementTypes
    # -------------------------------------------------------------------------
    def build_car_rep_map(body_style, label: str, is_suv: bool = False):
        length = 4.85 if is_suv else 4.65
        width = 1.94 if is_suv else 1.86
        lower_h = 0.72 if is_suv else 0.62
        cabin_l = 2.75 if is_suv else 2.35
        cabin_h = 0.64 if is_suv else 0.54
        wheel_z = 0.0
        items = [
            make_box_solid(length, width, lower_h, 0.0, 0.0, 0.28, body_style, f"{label} Chassis"),
            make_box_solid(cabin_l, width * 0.84, cabin_h, -0.22, 0.0, 0.28 + lower_h, style_cabin_glass, f"{label} Cabin"),
            make_box_solid(cabin_l * 0.82, width * 0.80, 0.08, -0.22, 0.0, 0.28 + lower_h + cabin_h, body_style, f"{label} Roof"),
            make_box_solid(0.72, width * 1.02, 0.64, length * 0.30, 0.0, wheel_z, style_tire_black, f"{label} Front Axle"),
            make_box_solid(0.72, width * 1.02, 0.64, -length * 0.30, 0.0, wheel_z, style_tire_black, f"{label} Rear Axle"),
        ]
        return make_rep_map(items), (length * width * 2 + 28.0, length * width * (lower_h + cabin_h * 0.5))

    def build_semi_truck_rep_map():
        items = [
            make_box_solid(6.2, 2.48, 1.35, 5.2, 0.0, 0.45, style_car_blue, "Tractor Chassis & Hood"),
            make_box_solid(2.8, 2.46, 1.95, 4.4, 0.0, 1.80, style_car_blue, "Tractor Sleeper Cab"),
            make_box_solid(0.4, 2.30, 0.65, 5.7, 0.0, 1.95, style_cabin_glass, "Tractor Windshield"),
            make_box_solid(13.6, 2.58, 2.85, -3.6, 0.0, 1.20, style_deck_white, "53-ft Reefer Trailer Body"),
            make_box_solid(13.62, 2.60, 0.35, -3.6, 0.0, 2.40, style_container_teal, "Trailer Side Stripe"),
            make_box_solid(0.95, 2.56, 0.95, 6.8, 0.0, 0.0, style_tire_black, "Steer Axle"),
            make_box_solid(2.20, 2.56, 0.95, 3.0, 0.0, 0.0, style_tire_black, "Drive Tandem Axles"),
            make_box_solid(2.20, 2.56, 0.95, -8.4, 0.0, 0.0, style_tire_black, "Trailer Tandem Axles"),
        ]
        return make_rep_map(items), (215.0, 118.5)

    def build_patrol_tow_rep_map():
        items = [
            make_box_solid(7.8, 2.42, 1.15, 0.0, 0.0, 0.38, style_safety_yellow, "CBBT Patrol Service Body"),
            make_box_solid(2.6, 2.32, 1.10, 1.6, 0.0, 1.53, style_safety_yellow, "CBBT Patrol Cab"),
            make_box_solid(1.1, 1.95, 0.35, 1.6, 0.0, 2.63, style_safety_orange, "CBBT Emergency Lightbar"),
            make_box_solid(3.2, 1.40, 1.15, -2.0, 0.0, 1.53, style_safety_orange, "Tunnel Tow Boom"),
            make_box_solid(0.88, 2.48, 0.85, 2.5, 0.0, 0.0, style_tire_black, "Front Axle"),
            make_box_solid(1.85, 2.48, 0.85, -2.2, 0.0, 0.0, style_tire_black, "Rear Dual Axle"),
        ]
        return make_rep_map(items), (86.4, 34.2)

    def build_container_ship_rep_map():
        hull_pts = [
            (-102.0, -15.0), (68.0, -16.0), (94.0, -10.0), (108.0, 0.0),
            (94.0, 10.0), (68.0, 16.0), (-102.0, 15.0), (-106.0, 10.0), (-106.0, -10.0)
        ]
        items = [
            make_poly_solid(hull_pts, 3.2, -1.0, style_antifouling, "Panamax Waterline Boot-Top Hull"),
            make_poly_solid(hull_pts, 9.5, 2.2, style_navy_hull, "Panamax Upper Freeboard Hull"),
            make_box_solid(62.0, 26.0, 7.8, 24.0, 0.0, 11.7, style_container_orange, "Forward Intermodal Container Bays"),
            make_box_solid(54.0, 26.0, 7.8, -42.0, 0.0, 11.7, style_container_teal, "Aft Intermodal Container Bays"),
            make_box_solid(48.0, 22.0, 5.2, 22.0, 0.0, 19.5, style_container_teal, "Upper Tier Forward Containers"),
            make_box_solid(18.0, 24.0, 16.5, -82.0, 0.0, 11.7, style_deck_white, "Navigation Bridge Superstructure"),
            make_box_solid(6.0, 28.0, 3.8, -80.0, 0.0, 28.2, style_cabin_glass, "Navigation Bridge Wings & Windows"),
            make_box_solid(7.5, 5.5, 8.0, -86.0, 0.0, 28.2, style_container_orange, "Exhaust Funnel Casing"),
        ]
        return make_rep_map(items), (18450.0, 94200.0)

    def build_naval_cargo_ship_rep_map():
        hull_pts = [
            (-84.0, -13.0), (54.0, -14.0), (78.0, -7.5), (90.0, 0.0),
            (78.0, 7.5), (54.0, 14.0), (-84.0, 13.0), (-88.0, 8.0), (-88.0, -8.0)
        ]
        items = [
            make_poly_solid(hull_pts, 2.8, -0.8, style_antifouling, "Naval Auxiliary Waterline Hull"),
            make_poly_solid(hull_pts, 8.8, 2.0, style_naval_grey, "Naval Auxiliary Main Hull"),
            make_box_solid(36.0, 21.0, 11.5, 22.0, 0.0, 10.8, style_naval_grey, "Forward Superstructure & Bridge"),
            make_box_solid(8.0, 23.0, 3.2, 32.0, 0.0, 22.3, style_cabin_glass, "Pilot House & Radar Bridge"),
            make_box_solid(4.5, 4.5, 11.0, 20.0, 0.0, 22.3, style_galv_steel, "Integrated Radar & TACAN Mast"),
            make_box_solid(46.0, 22.0, 5.5, -36.0, 0.0, 10.8, style_naval_grey, "Aft Cargo Hold & Helicopter Deck"),
        ]
        return make_rep_map(items), (13200.0, 58400.0)

    def build_tugboat_rep_map(length: float = 32.0, beam: float = 9.6):
        hull_pts = [
            (-length * 0.48, -beam * 0.46), (length * 0.28, -beam * 0.50),
            (length * 0.46, -beam * 0.24), (length * 0.52, 0.0),
            (length * 0.46, beam * 0.24), (length * 0.28, beam * 0.50),
            (-length * 0.48, beam * 0.46),
        ]
        items = [
            make_poly_solid(hull_pts, 3.2, -0.4, style_navy_hull, "Vessel Hull"),
            make_box_solid(length * 0.42, beam * 0.68, 3.4, length * 0.05, 0.0, 2.8, style_deck_white, "Deckhouse Cabin"),
            make_box_solid(length * 0.22, beam * 0.58, 2.6, length * 0.10, 0.0, 6.2, style_cabin_glass, "Wheelhouse"),
            make_box_solid(1.8, 1.8, 3.5, -length * 0.12, 0.0, 6.2, style_safety_orange, "Marine Radar Mast & Stack"),
        ]
        return make_rep_map(items), (640.0, 1120.0)

    def build_lcs_signal_rep_map():
        items = [
            make_box_solid(0.45, 1.15, 1.15, 0.0, 0.0, 0.0, style_tire_black, "LCS Matrix Signal Backplate"),
            make_box_solid(0.12, 0.85, 0.85, 0.24, 0.0, 0.15, style_signal_green, "LCS LED Green Arrow Aspect"),
        ]
        return make_rep_map(items), (4.2, 0.78)

    def build_ohvd_sensor_rep_map():
        # Compact 1.85m pedestal mounted on the +8.05m island coping outside the ramp railing
        items = [
            make_box_solid(0.28, 0.28, 1.45, 0.0, 0.0, 0.0, style_safety_yellow, "Overheight Detector Steel Pedestal Post"),
            make_box_solid(0.42, 0.42, 0.42, 0.0, 0.0, 1.45, style_safety_orange, "Dual-Beam Laser Transceiver Head"),
        ]
        return make_rep_map(items), (2.8, 0.24)

    def build_cctv_camera_rep_map():
        items = [
            make_box_solid(0.24, 0.24, 1.4, 0.0, 0.0, 0.0, style_galv_steel, "CCTV Mast Riser"),
            make_box_solid(0.55, 0.36, 0.38, 0.15, 0.0, 1.4, style_deck_white, "PTZ Environmental Camera Housing"),
        ]
        return make_rep_map(items), (2.4, 0.18)

    def build_fire_hydrant_cabinet_rep_map():
        items = [
            make_box_solid(0.95, 0.38, 1.65, 0.0, 0.0, 0.0, style_signal_red, "Tunnel Fire Standpipe & Hose Cabinet"),
            make_box_solid(0.75, 0.06, 0.90, 0.0, 0.18, 0.55, style_deck_white, "Fire Cabinet Inspection Glazing"),
        ]
        return make_rep_map(items), (5.2, 0.60)

    def build_sos_alarm_box_rep_map():
        items = [
            make_box_solid(0.65, 0.32, 1.40, 0.0, 0.0, 0.0, style_sos_blue, "Tunnel SOS Emergency Call & Pull-Box Enclosure"),
            make_box_solid(0.28, 0.28, 0.24, 0.0, 0.0, 1.40, style_safety_orange, "SOS Strobe Beacon"),
        ]
        return make_rep_map(items), (3.4, 0.31)

    def build_sump_deluge_skid_rep_map():
        items = [
            make_box_solid(3.6, 2.1, 0.30, 0.0, 0.0, 0.0, style_galv_steel, "Equipment Structural Steel Skid Base"),
            make_box_solid(1.5, 1.6, 1.85, -0.8, 0.0, 0.30, style_jet_blue, "Duplex Submersible Sump Pump & Manifold"),
            make_box_solid(1.4, 1.5, 1.95, 0.9, 0.0, 0.30, style_signal_red, "Portal Foam-Water Deluge Valve & Control Cabinet"),
        ]
        return make_rep_map(items), (32.5, 9.8)

    car_blue_rm, car_blue_q = build_car_rep_map(style_car_blue, "Atlantic Blue Passenger Sedan", False)
    car_red_rm, car_red_q = build_car_rep_map(style_car_red, "Crimson Crossover SUV", True)
    car_silver_rm, car_silver_q = build_car_rep_map(style_car_silver, "Titanium Silver Touring Sedan", False)
    car_white_rm, car_white_q = build_car_rep_map(style_car_white, "Pearl White Highway SUV", True)
    semi_rm, semi_q = build_semi_truck_rep_map()
    patrol_rm, patrol_q = build_patrol_tow_rep_map()
    container_ship_rm, container_ship_q = build_container_ship_rep_map()
    naval_ship_rm, naval_ship_q = build_naval_cargo_ship_rep_map()
    tug_rm, tug_q = build_tugboat_rep_map(32.0, 9.6)
    charter_rm, charter_q = build_tugboat_rep_map(22.0, 6.8)
    lcs_rm, lcs_q = build_lcs_signal_rep_map()
    ohvd_rm, ohvd_q = build_ohvd_sensor_rep_map()
    cctv_rm, cctv_q = build_cctv_camera_rep_map()
    fire_rm, fire_q = build_fire_hydrant_cabinet_rep_map()
    sos_rm, sos_q = build_sos_alarm_box_rep_map()
    sump_rm, sump_q = build_sump_deluge_skid_rep_map()

    def create_type(ifc_type_cls: str, name: str, desc: str, rep_map, predef: str):
        return f.create_entity(
            ifc_type_cls,
            GlobalId=new_guid(),
            OwnerHistory=owner_history,
            Name=name,
            Description=desc,
            RepresentationMaps=[rep_map],
            PredefinedType=predef,
        )

    type_car_blue = create_type("IfcVehicleType", "Passenger Sedan Type A (Atlantic Metallic Blue)", "4-door highway passenger vehicle on US Route 13 CBBT corridor", car_blue_rm, "VEHICLEWHEELED")
    type_car_red = create_type("IfcVehicleType", "Crossover SUV Type B (Crimson Red)", "All-wheel-drive sport utility vehicle on US Route 13 CBBT corridor", car_red_rm, "VEHICLEWHEELED")
    type_car_silver = create_type("IfcVehicleType", "Touring Sedan Type C (Titanium Silver)", "4-door touring sedan on US Route 13 CBBT corridor", car_silver_rm, "VEHICLEWHEELED")
    type_car_white = create_type("IfcVehicleType", "Highway SUV Type D (Pearl White)", "Full-size sport utility vehicle on US Route 13 CBBT corridor", car_white_rm, "VEHICLEWHEELED")
    type_semi = create_type("IfcVehicleType", "WB-67 18-Wheeler Tractor-Trailer Freight Unit Type", "AASHTO WB-67 5-axle interstate tractor-semitrailer freight vehicle", semi_rm, "VEHICLEWHEELED")
    type_patrol = create_type("IfcVehicleType", "CBBT Authority Emergency Patrol & Tunnel Recovery Tow Unit Type", "CBBT Police & Tunnel Operations heavy-duty emergency response & recovery vehicle", patrol_rm, "VEHICLEWHEELED")

    type_container_ship = create_type("IfcVehicleType", "Panamax Container Vessel Type (210m LOA)", "Deep-draft ocean-going container ship transiting Thimble Shoal Channel above immersed tunnel", container_ship_rm, "VEHICLEMARINE")
    type_naval_ship = create_type("IfcVehicleType", "Naval Auxiliary Cargo Vessel Type (175m LOA)", "Deep-draft naval auxiliary & logistics vessel transiting Chesapeake Channel above immersed tunnel", naval_ship_rm, "VEHICLEMARINE")
    type_tug = create_type("IfcVehicleType", "Harbor & Channel Escort Tugboat Type (32m LOA)", "Twin-screw Z-drive marine escort tugboat in North Channel span", tug_rm, "VEHICLEMARINE")
    type_charter = create_type("IfcVehicleType", "Chesapeake Bay Pilot & Charter Boat Type (22m LOA)", "Coastal pilot & sport-fishing vessel at Sea Gull Pier", charter_rm, "VEHICLEMARINE")

    type_lcs = create_type("IfcSignalType", "CBBT Overhead LED Lane-Control Signal (LCS) Matrix Type", "MUTCD 4M reversible/tunnel lane-control signal head (Green Downward Arrow / Red X)", lcs_rm, "VISUAL")
    type_ohvd = create_type("IfcSensorType", "Dual-Beam Laser Overheight Vehicle Detector (OHVD) Type", "Dual photoelectric/laser overheight vehicle detector protecting 13ft-6in tunnel portals", ohvd_rm, "LIGHTSENSOR")
    type_cctv = create_type("IfcAudioVisualApplianceType", "Marine-Grade PTZ Traffic & Tunnel Surveillance Camera Type", "NEMA 4X IP68 pan-tilt-zoom HD traffic & incident detection camera", cctv_rm, "CAMERA")
    type_fire = create_type("IfcFireSuppressionTerminalType", "Tunnel Standpipe & Fire Hose Cabinet Niche Assembly Type", "NFPA 502 Class I/III dual 2.5-in fire department standpipe valve and hose cabinet", fire_rm, "FIREHYDRANT")
    type_sos = create_type("IfcAlarmType", "Tunnel SOS Emergency Telephone & Fire Pull-Box Station Type", "NFPA 502 blue-light emergency motorist call box and manual fire alarm pull station", sos_rm, "MANUALPULLBOX")
    type_sump = create_type("IfcUnitaryEquipmentType", "Tunnel Portal Sump Drainage & Deluge Valve Skid Type", "Duplex submersible portal drainage pump and foam-water deluge control skid", sump_rm, "AIRHANDLER")

    type_instances = {t.id(): [] for t in [
        type_car_blue, type_car_red, type_car_silver, type_car_white, type_semi, type_patrol,
        type_container_ship, type_naval_ship, type_tug, type_charter,
        type_lcs, type_ohvd, type_cctv, type_fire, type_sos, type_sump
    ]}

    # Map container ID -> IfcRelContainedInSpatialStructure and world origin
    cont_rels = {}
    cont_origins = {}
    cont_by_name = {}
    for r in f.by_type("IfcRelContainedInSpatialStructure"):
        cid = r.RelatingStructure.id()
        cname = r.RelatingStructure.Name or ""
        cont_rels[cid] = r
        cont_by_name[cname] = cid
        m = pl.get_local_placement(r.RelatingStructure.ObjectPlacement)
        cont_origins[cid] = (float(m[0, 3]), float(m[1, 3]), float(m[2, 3]), r.RelatingStructure.ObjectPlacement)

    def find_container(substr: str, default_cid: int = None) -> int:
        for cname, cid in cont_by_name.items():
            if substr in cname:
                return cid
        return default_cid if default_cid is not None else next(iter(cont_rels.keys()))

    cid_trestle_a = find_container("Trestles A/B/C")
    cid_south_plaza = find_container("Virginia Beach South Toll Plaza", cid_trestle_a)
    cid_north_plaza = find_container("Wise Point North Toll Plaza", cid_trestle_a)
    cid_isl1 = find_container("Island #1", cid_trestle_a)
    cid_isl2 = find_container("Island #2", cid_trestle_a)
    cid_isl3 = find_container("Island #3", cid_trestle_a)
    cid_isl4 = find_container("Island #4", cid_trestle_a)
    cid_pier = find_container("Sea Gull Fishing Pier", cid_isl1)
    cid_nc_sub = find_container("North Channel Substructure", cid_trestle_a)
    cid_dolphins = find_container("Shipping Channel Protection", cid_isl1)
    cid_thimble_arch = find_container("Thimble Shoal Immersed Tunnel — Structural", cid_isl1)
    cid_thimble_mep = find_container("Thimble Shoal Tunnel — Jet Fans", cid_thimble_arch)
    cid_ches_arch = find_container("Chesapeake Channel Immersed Tunnel — Structural", cid_isl3)
    cid_ches_mep = find_container("Chesapeake Channel Tunnel — Jet Fans", cid_ches_arch)

    new_steel_elements = []
    new_enclosure_elements = []

    def spawn_element(
        ifc_cls: str,
        el_type,
        rep_map,
        q_vals,
        container_id: int,
        name: str,
        desc: str,
        tag: str,
        predef: str,
        world_xyz,
        x_dir,
        z_up=(0.0, 0.0, 1.0),
        is_enclosure: bool = False,
    ):
        cx0, cy0, cz0, parent_lp = cont_origins[container_id]
        lx = float(world_xyz[0]) - cx0
        ly = float(world_xyz[1]) - cy0
        lz = float(world_xyz[2]) - cz0

        zv = np.array(z_up, dtype=float)
        zv /= np.linalg.norm(zv)
        xv = np.array(x_dir, dtype=float)
        xv = xv - np.dot(xv, zv) * zv
        if np.linalg.norm(xv) < 1e-6:
            xv = np.array([1.0, 0.0, 0.0])
        else:
            xv /= np.linalg.norm(xv)

        ap = f.createIfcAxis2Placement3D(
            f.createIfcCartesianPoint((round(lx, 4), round(ly, 4), round(lz, 4))),
            f.createIfcDirection((round(float(zv[0]), 6), round(float(zv[1]), 6), round(float(zv[2]), 6))),
            f.createIfcDirection((round(float(xv[0]), 6), round(float(xv[1]), 6), round(float(xv[2]), 6))),
        )
        lp = f.createIfcLocalPlacement(parent_lp, ap)
        mapped_item = f.createIfcMappedItem(rep_map, ident_xform)
        shape_rep = f.createIfcShapeRepresentation(body_ctx, "Body", "MappedRepresentation", [mapped_item])
        prod_shape = f.createIfcProductDefinitionShape(None, None, [shape_rep])

        el = f.create_entity(
            ifc_cls,
            GlobalId=new_guid(),
            OwnerHistory=owner_history,
            Name=name,
            Description=desc,
            ObjectType=el_type.Name,
            ObjectPlacement=lp,
            Representation=prod_shape,
            Tag=tag,
            PredefinedType=predef,
        )

        area, vol = q_vals
        q_area1 = f.createIfcQuantityArea("GrossSurfaceArea", None, None, round(float(area), 3), None)
        q_area2 = f.createIfcQuantityArea("NetSurfaceArea", None, None, round(float(area), 3), None)
        q_vol1 = f.createIfcQuantityVolume("GrossVolume", None, None, round(float(vol), 3), None)
        q_vol2 = f.createIfcQuantityVolume("NetVolume", None, None, round(float(vol), 3), None)
        qto = f.createIfcElementQuantity(new_guid(), owner_history, "Qto_BodyGeometryValidation", "Calculated physical geometry quantities", "IfcOpenShell", [q_area1, q_area2, q_vol1, q_vol2])
        f.createIfcRelDefinesByProperties(new_guid(), owner_history, None, None, [el], qto)

        p1 = f.createIfcPropertySingleValue("AgencyBridgeId", None, f.createIfcIdentifier("VA-CBBT-US13-CORRIDOR"), None)
        p2 = f.createIfcPropertySingleValue("NationalBridgeId", None, f.createIfcIdentifier("NBI-5100130CBBT01"), None)
        p3 = f.createIfcPropertySingleValue("BridgeGroup", None, f.createIfcLabel("CBBT Operations, ITS & Corridor Traffic Division"), None)
        p4 = f.createIfcPropertySingleValue("HistoricalSignificance", None, f.createIfcLabel("ASCE Seven Engineering Wonders of the Modern World — CBBT"), None)
        pset = f.createIfcPropertySet(new_guid(), owner_history, "usBrPset_BridgeIdentification", None, [p1, p2, p3, p4])
        f.createIfcRelDefinesByProperties(new_guid(), owner_history, None, None, [el], pset)

        type_instances[el_type.id()].append(el)
        cont_rels[container_id].RelatedElements = list(cont_rels[container_id].RelatedElements) + [el]
        if is_enclosure:
            new_enclosure_elements.append(el)
        else:
            new_steel_elements.append(el)
        return el

    # -------------------------------------------------------------------------
    # 2. Place Wheeled Vehicles along Mainline IfcPavement & Tunnel Roadways
    # -------------------------------------------------------------------------
    elem_to_container = {}
    for cid, rel in cont_rels.items():
        for e in rel.RelatedElements:
            elem_to_container[e.id()] = cid

    pave_by_tag = {}
    oriented_road_segs = []
    tunnel_road_segs = []
    for e in f.by_type("IfcPavement"):
        tag = getattr(e, "Tag", "") or ""
        m = pl.get_local_placement(e.ObjectPlacement)
        it = e.Representation.Representations[0].Items[0]
        depth = float(getattr(it, "Depth", 0.28))
        prof = getattr(it, "SweptArea", None)
        seg_len = float(getattr(prof, "XDim", 20.0)) if prof else 20.0
        cid = elem_to_container.get(e.id(), cid_trestle_a)
        pave_by_tag[tag] = (e.Name or "Roadway", m, depth, seg_len, cid)
        if tag.startswith("TPAVE-"):
            tunnel_road_segs.append((e.Name or "Tunnel Roadway", m, depth, cid, tag))
        elif tag.startswith("PAVE-"):
            # Only include primary carriageways & open-cut ramps (exclude short transitions < 12m)
            if seg_len >= 12.0:
                oriented_road_segs.append((e.Name or "Roadway", m, depth, seg_len, cid, tag))

    # Also check IfcBuiltElement in case TPAVE is IfcBuiltElement
    if not tunnel_road_segs:
        for e in f.by_type("IfcBuiltElement"):
            if "Internal Asphalt Roadway Wearing Surface" in (e.Name or ""):
                tag = getattr(e, "Tag", "") or ""
                m = pl.get_local_placement(e.ObjectPlacement)
                it = e.Representation.Representations[0].Items[0]
                depth = float(getattr(it, "Depth", 0.28))
                cid = elem_to_container.get(e.id(), cid_thimble_arch)
                tunnel_road_segs.append((e.Name, m, depth, cid, tag))

    car_types = [
        (type_car_blue, car_blue_rm, car_blue_q, "Atlantic Blue Passenger Sedan"),
        (type_car_red, car_red_rm, car_red_q, "Crimson Crossover SUV"),
        (type_car_silver, car_silver_rm, car_silver_q, "Titanium Silver Touring Sedan"),
        (type_car_white, car_white_rm, car_white_q, "Pearl White Highway SUV"),
    ]

    veh_count = 0
    used_pave_tags = set()

    # First, explicitly select 8 Portal Island surface roadway & sloped open-cut ramp segments
    # so every island (#1, #2, #3, #4) has vehicles both on its flat surface roadway AND on its sloped open-cut ramp,
    # all computed directly from the segment's 3D placement matrix `m` (zero floating)!
    island_featured_tags = [
        "PAVE-NB-90583844-2",  # Island #1 Surface Transition Roadway (z = 8.50m)
        "PAVE-NB-41336138-3",  # Island #1 Open-Cut Depressed Ramp (sloped, z = 4.19m)
        "PAVE-NB-90583872-3",  # Island #2 Open-Cut Depressed Ramp (sloped, z = 2.71m)
        "PAVE-NB-20375737-2",  # Island #2 Surface Transition Roadway (z = 8.50m)
        "PAVE-NB-90583865-2",  # Island #3 Surface Transition Roadway (z = 8.50m)
        "PAVE-NB-46930465-2",  # Island #3 Open-Cut Depressed Ramp (sloped, z = 5.57m)
        "PAVE-NB-90583871-3",  # Island #4 Open-Cut Depressed Ramp (sloped, z = 2.40m)
        "PAVE-NB-46930462-2",  # Island #4 Surface Transition Roadway (z = 8.50m)
    ]

    # Sample 48 additional corridor segments across the 28.6 km bridge-tunnel
    remaining_segs = [s for s in oriented_road_segs if s[5] not in island_featured_tags]
    step = max(1, len(remaining_segs) // 48)
    sampled_corridor_segs = remaining_segs[::step][:48]
    for ft in island_featured_tags:
        if ft in pave_by_tag:
            sname, m, depth, seg_len, cid = pave_by_tag[ft]
            sampled_corridor_segs.append((sname, m, depth, seg_len, cid, ft))

    for idx, (seg_name, m, depth, seg_len, cid, tag) in enumerate(sampled_corridor_segs):
        used_pave_tags.add(tag)
        is_sb = "(SB)" in seg_name
        tangent = m[:3, 0].copy()
        left_vec = m[:3, 1].copy()
        up_vec = m[:3, 2].copy()
        travel_dir = -tangent if is_sb else tangent
        lane_sign = 1.0 if (idx % 2 == 0) else -1.0
        lateral_offset = lane_sign * 2.05
        pos = m[:3, 3] + up_vec * depth + left_vec * lateral_offset

        veh_count += 1
        direction_tag = "SB" if is_sb else "NB"
        # Only place 18-Wheeler semi-trucks on long straight/gentle segments (seg_len >= 26.0m, excluding steep open-cut ramps)
        if (idx % 5 == 2) and (seg_len >= 26.0) and (tag not in island_featured_tags) and ("Open-Cut" not in seg_name):
            spawn_element(
                "IfcVehicle", type_semi, semi_rm, semi_q, cid,
                f"US-13 {direction_tag} Interstate Freight Tractor-Trailer #{veh_count:02d} ({seg_name.split('(')[0].strip()})",
                f"5-axle WB-67 tractor-semitrailer traveling {direction_tag} on {seg_name}",
                f"CBBT-VEH-{veh_count:03d}", "VEHICLEWHEELED",
                pos, travel_dir, up_vec, is_enclosure=False
            )
        else:
            ct, crm, cq, cname = car_types[idx % len(car_types)]
            spawn_element(
                "IfcVehicle", ct, crm, cq, cid,
                f"US-13 {direction_tag} {cname} #{veh_count:02d} ({seg_name.split('(')[0].strip()})",
                f"Passenger vehicle traveling {direction_tag} along {seg_name}",
                f"CBBT-VEH-{veh_count:03d}", "VEHICLEWHEELED",
                pos, travel_dir, up_vec, is_enclosure=False
            )

    # Place subaqueous tunnel interior vehicles in the deep channel segments (excluding portal-edge segments under the islands)
    deep_tunnel_segs = [
        s for s in tunnel_road_segs
        if not s[4].endswith("-01") and not s[4].endswith("-02")
        and not s[4].endswith("-13") and not s[4].endswith("-14")
    ]
    for tidx, (tname, m, depth, cid, tag) in enumerate(deep_tunnel_segs[::3]):
        tangent = m[:3, 0].copy()
        left_vec = m[:3, 1].copy()
        up_vec = m[:3, 2].copy()
        is_sb_lane = (tidx % 2 == 1)
        travel_dir = -tangent if is_sb_lane else tangent
        lane_offset = 1.95 if is_sb_lane else -1.95
        pos = m[:3, 3] + up_vec * depth + left_vec * lane_offset
        veh_count += 1
        ct, crm, cq, cname = car_types[tidx % len(car_types)]
        spawn_element(
            "IfcVehicle", ct, crm, cq, cid,
            f"Tunnel Interior {'SB' if is_sb_lane else 'NB'} Lane — {cname} #{veh_count:02d} ({tname.split('1964')[0].strip()} Tunnel)",
            f"Motor vehicle transiting submerged tube section of {tname}",
            f"CBBT-VEH-{veh_count:03d}", "VEHICLEWHEELED",
            pos, travel_dir, up_vec, is_enclosure=False
        )

    # Place 6 CBBT Authority Emergency Patrol & Recovery Tow Trucks directly on the paved Island Pull-Off Shoulder Aprons & Toll Plazas
    # Using exact pavement placement matrices `m` + `up_vec * depth` and aligned parallel to the adjacent roadway!
    patrol_apron_specs = [
        (cid_south_plaza, "Virginia Beach South Toll Plaza — CBBT Police & Emergency Patrol Unit #1", "PAVE-NB-288016838-2", None, -2.2, 1.0),
        (cid_isl1, "Island #1 (South Thimble Shoal Portal) — Tunnel Emergency Recovery Tow Truck #1", "PPAVE-NB-90583844-3", (0.2896, 0.9571, 0.0), 0.0, 1.0),
        (cid_isl2, "Island #2 (North Thimble Shoal Portal) — Tunnel Emergency Recovery Tow Truck #2", "PPAVE-NB-20375737-3", (-0.2896, -0.9571, 0.0), 0.0, 1.0),
        (cid_isl3, "Island #3 (South Chesapeake Channel Portal) — Tunnel Emergency Recovery Tow Truck #3", "PPAVE-NB-90583865-3", (0.7406, 0.6720, 0.0), 0.0, 1.0),
        (cid_isl4, "Island #4 (North Chesapeake Channel Portal) — Tunnel Emergency Recovery Tow Truck #4", "PPAVE-NB-46930462-3", (-0.7406, -0.6720, 0.0), 0.0, 1.0),
        (cid_north_plaza, "Eastern Shore Wise Point North Toll Plaza — CBBT Police & Emergency Patrol Unit #2", "PAVE-SB-20375685-2", None, -2.2, -1.0),
    ]
    for cid, pname, ptag, explicit_dir, lat_off, dir_sign in patrol_apron_specs:
        veh_count += 1
        if ptag in pave_by_tag:
            _, pm, pdepth, _, _ = pave_by_tag[ptag]
            if explicit_dir is not None:
                ptangent = np.array(explicit_dir, dtype=float)
            else:
                ptangent = pm[:3, 0].copy() * float(dir_sign)
            pleft = pm[:3, 1].copy()
            pup = pm[:3, 2].copy()
            pxyz = pm[:3, 3] + pup * pdepth + pleft * float(lat_off)
            spawn_element(
                "IfcVehicle", type_patrol, patrol_rm, patrol_q, cid,
                pname,
                f"CBBT District Tunnel Operations standby emergency response and motorist assistance tow vehicle parked on {ptag}",
                f"CBBT-VEH-{veh_count:03d}", "VEHICLEWHEELED",
                pxyz, ptangent, pup, is_enclosure=False
            )

    # -------------------------------------------------------------------------
    # 3. Place Marine Vessels (IfcVehicle, VEHICLEMARINE) in the Shipping Channels
    # -------------------------------------------------------------------------
    marine_vessels = [
        (
            cid_dolphins, type_container_ship, container_ship_rm, container_ship_q,
            "Thimble Shoal Federal Navigation Channel — Panamax Container Ship 'MV Chesapeake Bay Express'",
            "210m Panamax container vessel outbound from Port of Virginia crossing directly over the Thimble Shoal Immersed-Tube Tunnel",
            "CBBT-SHIP-001",
            (1791.8, 6848.2, 0.0),
            (0.9571, -0.2896, 0.0),
        ),
        (
            cid_dolphins, type_naval_ship, naval_ship_rm, naval_ship_q,
            "Chesapeake Main Navigation Channel — Naval Auxiliary Vessel 'USNS Hampton Roads'",
            "175m Military Sealift Command auxiliary cargo vessel crossing above the Chesapeake Channel Immersed-Tube Tunnel",
            "CBBT-SHIP-002",
            (5132.6, 14417.0, 0.0),
            (-0.6720, 0.7406, 0.0),
        ),
        (
            cid_nc_sub, type_tug, tug_rm, tug_q,
            "North Channel High-Level Bridge Navigation Span — Marine Escort Tugboat 'MV Cape Charles'",
            "32m twin-screw harbor and channel escort tugboat transiting beneath the North Channel High-Level Bridge main truss span",
            "CBBT-SHIP-003",
            (12284.8, 19947.3, 0.0),
            (-0.3649, 0.9311, 0.0),
        ),
        (
            cid_pier, type_charter, charter_rm, charter_q,
            "Sea Gull Fishing Pier (Island #1) — Chesapeake Bay Pilot & Charter Vessel 'MV Sea Gull II'",
            "22m coastal pilot and sport-fishing vessel moored alongside the Sea Gull Pier marine berth at South Thimble Island",
            "CBBT-SHIP-004",
            (1328.0, 6155.0, 0.0),
            (0.2896, 0.9571, 0.0),
        ),
    ]
    for cid, stype, srm, sq, sname, sdesc, stag, sxyz, sdir in marine_vessels:
        spawn_element(
            "IfcVehicle", stype, srm, sq, cid,
            sname, sdesc, stag, "VEHICLEMARINE",
            sxyz, sdir, (0.0, 0.0, 1.0), is_enclosure=False
        )

    # -------------------------------------------------------------------------
    # 4. Place Tunnel Portal Lane-Control Signals (IfcSignal), Overheight Laser
    #    Detectors (IfcSensor), and PTZ CCTV Cameras (IfcAudioVisualAppliance)
    # -------------------------------------------------------------------------
    portal_configs = [
        (cid_thimble_mep, "Island #1 (South Thimble Shoal Portal)", (1496.63, 6028.08, 8.05), (0.2896, 0.9571, 0.0), (-0.9571, 0.2896, 0.0)),
        (cid_thimble_mep, "Island #2 (North Thimble Shoal Portal)", (2003.36, 7702.73, 8.05), (-0.2896, -0.9571, 0.0), (0.9571, -0.2896, 0.0)),
        (cid_ches_mep, "Island #3 (South Chesapeake Channel Portal)", (4492.65, 13897.05, 8.05), (0.7406, 0.6720, 0.0), (-0.6720, 0.7406, 0.0)),
        (cid_ches_mep, "Island #4 (North Chesapeake Channel Portal)", (5722.06, 15012.57, 8.05), (-0.7406, -0.6720, 0.0), (0.6720, -0.7406, 0.0)),
    ]

    lcs_idx = 0
    ohvd_idx = 0
    cctv_idx = 0
    sump_idx = 0

    for cid, pname, pxyz, fdir, ldir in portal_configs:
        px, py, pz_deck = pxyz
        fx, fy, _ = fdir
        lx, ly, _ = ldir
        approach_facing = (-fx, -fy, 0.0)

        # Two Lane-Control Signals (NB & SB lanes) mounted on RampCrossStrut_2 (42m ahead of portal at z = 6.95m) facing approaching traffic
        for lane_label, lat_off in [("Northbound Lane", -1.95), ("Southbound Lane", 1.95)]:
            lcs_idx += 1
            spawn_element(
                "IfcSignal", type_lcs, lcs_rm, lcs_q, cid,
                f"{pname} — Approach Ramp {lane_label} LED Lane-Control Signal #{lcs_idx:02d}",
                f"MUTCD fiber-optic/LED lane-control signal head mounted on overhead ramp strut at {pname}",
                f"CBBT-LCS-{lcs_idx:03d}", "VISUAL",
                (px - fx * 42.0 + lx * lat_off, py - fy * 42.0 + ly * lat_off, 6.85),
                approach_facing, (0.0, 0.0, 1.0), is_enclosure=True
            )

        # Dual Overheight Vehicle Laser Detectors seated directly on the +8.05m Island Concrete Deck outside the ramp railings (42m ahead of portal)
        for side_label, lat_off in [("East Roadside Transceiver", -5.2), ("West Roadside Receiver", 5.2)]:
            ohvd_idx += 1
            spawn_element(
                "IfcSensor", type_ohvd, ohvd_rm, ohvd_q, cid,
                f"{pname} Approach — 13ft-6in Overheight Vehicle Laser Detector ({side_label}) #{ohvd_idx:02d}",
                f"Dual-beam infrared/laser overheight vehicle detection post seated on +8.05m island deck at {pname}",
                f"CBBT-OHVD-{ohvd_idx:03d}", "LIGHTSENSOR",
                (px - fx * 42.0 + lx * lat_off, py - fy * 42.0 + ly * lat_off, pz_deck),
                fdir,
                (0.0, 0.0, 1.0), is_enclosure=True
            )

        # Portal High-Definition PTZ Traffic & Marine Surveillance Camera mounted on the +12.40m Ventilation Building Parapet Roof
        cctv_idx += 1
        spawn_element(
            "IfcAudioVisualAppliance", type_cctv, cctv_rm, cctv_q, cid,
            f"{pname} — High-Definition PTZ Tunnel Approach & Channel Surveillance Camera #{cctv_idx:02d}",
            f"Pan-tilt-zoom IP68 marine traffic surveillance camera on ventilation building roof at {pname}",
            f"CBBT-CCTV-{cctv_idx:03d}", "CAMERA",
            (px - fx * 2.5 + lx * 5.2, py - fy * 2.5 + ly * 5.2, 12.40),
            approach_facing, (0.0, 0.0, 1.0), is_enclosure=True
        )

        # Portal Sump Drainage & Foam-Water Deluge Control Skid seated on the +8.05m Island Concrete Deck beside the ventilation building
        sump_idx += 1
        spawn_element(
            "IfcUnitaryEquipment", type_sump, sump_rm, sump_q, cid,
            f"{pname} — Portal Sump Drainage Pump & Deluge Fire Control Skid #{sump_idx:02d}",
            f"Duplex storm/seawater portal sump pump and deluge valve control skid seated on +8.05m island deck at {pname}",
            f"CBBT-SUMP-{sump_idx:03d}", "AIRHANDLER",
            (px + fx * 12.0 - lx * 10.5, py + fy * 12.0 - ly * 10.5, pz_deck),
            fdir, (0.0, 0.0, 1.0), is_enclosure=False
        )

    # Also mount Lane-Control Signals & PTZ Cameras directly on the actual Overhead Sign Gantry Crossbeams (`GBEAM-01`, `GBEAM-02`, `GBEAM-16`)
    gbeam_by_tag = {getattr(b, "Tag", ""): b for b in f.by_type("IfcBeam") if (getattr(b, "Tag", "") or "").startswith("GBEAM-")}
    for g_tag in ("GBEAM-01", "GBEAM-02", "GBEAM-16"):
        if g_tag in gbeam_by_tag:
            gb = gbeam_by_tag[g_tag]
            gm = pl.get_local_placement(gb.ObjectPlacement)
            gname = (gb.Name or g_tag).split("—")[0].strip()
            gcid = elem_to_container.get(gb.id(), cid_south_plaza)
            gx, gy, gz = float(gm[0, 3]), float(gm[1, 3]), float(gm[2, 3])
            gtangent = gm[:3, 0].copy()
            gleft = gm[:3, 1].copy()
            is_sb_gantry = "SB" in gname
            g_facing = gtangent if is_sb_gantry else -gtangent
            for lane_no, lat_off in [(1, -2.05), (2, 2.05)]:
                lcs_idx += 1
                spawn_element(
                    "IfcSignal", type_lcs, lcs_rm, lcs_q, gcid,
                    f"{gname} — Lane {lane_no} Overhead LED Lane-Control Signal #{lcs_idx:02d}",
                    f"Overhead gantry LED lane-control signal aspect mounted on {gname}",
                    f"CBBT-LCS-{lcs_idx:03d}", "VISUAL",
                    (gx + float(gleft[0]) * lat_off, gy + float(gleft[1]) * lat_off, gz - 0.35),
                    g_facing, (0.0, 0.0, 1.0), is_enclosure=True
                )
            cctv_idx += 1
            spawn_element(
                "IfcAudioVisualAppliance", type_cctv, cctv_rm, cctv_q, gcid,
                f"{gname} — Overhead PTZ Traffic Monitoring Camera #{cctv_idx:02d}",
                f"Overhead gantry pan-tilt-zoom traffic camera mounted on {gname}",
                f"CBBT-CCTV-{cctv_idx:03d}", "CAMERA",
                (gx, gy, gz + 0.78),
                g_facing, (0.0, 0.0, 1.0), is_enclosure=True
            )

    # -------------------------------------------------------------------------
    # 5. Place Tunnel Interior Fire Standpipe Cabinets (IfcFireSuppressionTerminal),
    #    Blue SOS Emergency Call Boxes (IfcAlarm), and Tunnel Bore CCTV Cameras
    # -------------------------------------------------------------------------
    fire_idx = 0
    sos_idx = 0
    for tidx, (tname, m, depth, cid, tag) in enumerate(tunnel_road_segs[::2]):
        tangent = m[:3, 0].copy()
        left_vec = m[:3, 1].copy()
        up_vec = m[:3, 2].copy()
        mep_cid = cid_thimble_mep if cid == cid_thimble_arch else cid_ches_mep
        short_tunnel = tname.split("Internal")[0].strip()

        fire_pos = m[:3, 3] + up_vec * (depth + 0.35) - left_vec * 4.15
        fire_idx += 1
        spawn_element(
            "IfcFireSuppressionTerminal", type_fire, fire_rm, fire_q, mep_cid,
            f"{short_tunnel} — Fire Standpipe & Hose Cabinet Niche #{fire_idx:02d}",
            f"NFPA 502 Class I/III fire department standpipe valve and 150-ft hose reel cabinet in {short_tunnel}",
            f"CBBT-FIRE-{fire_idx:03d}", "FIREHYDRANT",
            fire_pos, tangent, up_vec, is_enclosure=True
        )

        sos_pos = m[:3, 3] + up_vec * (depth + 0.35) + left_vec * 4.15
        sos_idx += 1
        spawn_element(
            "IfcAlarm", type_sos, sos_rm, sos_q, mep_cid,
            f"{short_tunnel} — Blue SOS Emergency Telephone & Fire Alarm Pull-Box #{sos_idx:02d}",
            f"Direct-line tunnel control room SOS emergency telephone and manual fire alarm pull station in {short_tunnel}",
            f"CBBT-SOS-{sos_idx:03d}", "MANUALPULLBOX",
            sos_pos, -tangent, up_vec, is_enclosure=True
        )

        if tidx % 2 == 0:
            cctv_pos = m[:3, 3] + up_vec * 3.85 - left_vec * 3.80
            cctv_idx += 1
            spawn_element(
                "IfcAudioVisualAppliance", type_cctv, cctv_rm, cctv_q, mep_cid,
                f"{short_tunnel} — Interior Bore Incident Detection CCTV Camera #{cctv_idx:02d}",
                f"Automatic incident detection (AID) & thermal tunnel bore camera in {short_tunnel}",
                f"CBBT-CCTV-{cctv_idx:03d}", "CAMERA",
                cctv_pos, tangent, up_vec, is_enclosure=True
            )

    # -------------------------------------------------------------------------
    # 6. Register IfcRelDefinesByType, IfcRelAssociatesMaterial & Classifications
    # -------------------------------------------------------------------------
    for tid, els in type_instances.items():
        if els:
            t_ent = f.by_id(tid)
            f.createIfcRelDefinesByType(new_guid(), owner_history, f"Type Assignment for {t_ent.Name}", None, els, t_ent)

    all_new_elements = new_steel_elements + new_enclosure_elements

    for mrel in f.by_type("IfcRelAssociatesMaterial"):
        mname = getattr(mrel.RelatingMaterial, "Name", "") or ""
        if "Weathering" in mname or "StructuralSteel" in mname:
            mrel.RelatedObjects = list(mrel.RelatedObjects) + new_steel_elements
            break

    for mrel in f.by_type("IfcRelAssociatesMaterial"):
        mname = getattr(mrel.RelatingMaterial, "Name", "") or ""
        if "Louver" in mname or "Aluminum" in mname:
            mrel.RelatedObjects = list(mrel.RelatedObjects) + new_enclosure_elements
            break

    for crel in f.by_type("IfcRelAssociatesClassification"):
        c_ident = getattr(crel.RelatingClassification, "Identification", "") or ""
        if c_ident == "L-PHD":
            crel.RelatedObjects = list(crel.RelatedObjects) + all_new_elements
        elif c_ident == "usBridge_Utilities":
            crel.RelatedObjects = list(crel.RelatedObjects) + all_new_elements
        elif c_ident == "RIH-STR-STEEL":
            crel.RelatedObjects = list(crel.RelatedObjects) + new_steel_elements
        elif c_ident == "RIH-ARCH-ENC":
            crel.RelatedObjects = list(crel.RelatedObjects) + new_enclosure_elements

    f.write(IFC_PATH)
    print(f"Successfully added {len(all_new_elements)} new IFC4X3 elements to {IFC_PATH}:")
    print(f"  - {veh_count} Wheeled Vehicles (IfcVehicle, VEHICLEWHEELED: Cars, SUVs, 18-Wheelers, CBBT Patrol Tow Trucks)")
    print(f"  - {len(marine_vessels)} Marine Vessels (IfcVehicle, VEHICLEMARINE: Container Ship, Naval Cargo Ship, Escort Tug, Charter Boat)")
    print(f"  - {lcs_idx} Overhead LED Lane-Control Signals (IfcSignal, VISUAL)")
    print(f"  - {ohvd_idx} Overheight Vehicle Laser Detectors (IfcSensor, LIGHTSENSOR)")
    print(f"  - {cctv_idx} PTZ Traffic & Tunnel Surveillance Cameras (IfcAudioVisualAppliance, CAMERA)")
    print(f"  - {fire_idx} Tunnel Fire Standpipe & Hose Cabinets (IfcFireSuppressionTerminal, FIREHYDRANT)")
    print(f"  - {sos_idx} Tunnel Blue SOS Emergency Call Stations (IfcAlarm, MANUALPULLBOX)")
    print(f"  - {sump_idx} Portal Sump & Deluge Control Skids (IfcUnitaryEquipment, AIRHANDLER)")


if __name__ == "__main__":
    main()
