# Chesapeake Bay Bridge-Tunnel — IFC 4.3 test model

`chesapeake_bay_bridge_tunnel_4x3.ifc` is a large, multi-facility IFC 4.3 model of the Chesapeake Bay Bridge-Tunnel (CBBT) corridor in Virginia, USA. It was generated with parametric Python scripts (IfcOpenShell) using open geodata and photo references. Its purpose is to test and demonstrate forensicBIM audits on infrastructure models.

> **This is a synthetic model, not an as-built record.** It was not produced, reviewed or endorsed by the Chesapeake Bay Bridge and Tunnel District, VDOT or FHWA. Geometry is approximate. The asset, inspection, warranty and design-loading properties are placeholder values. Do not use it for engineering, inspection or maintenance decisions.

## File summary

| Item | Value |
| :--- | :--- |
| Schema | `IFC4X3_ADD2` |
| Model view | ReferenceView, Alignment-basedView |
| Authoring | IfcOpenShell 0.8.0 ("Forensic-BIM Civil & Marine Infrastructure Authoring Engine 2.0") |
| File size | ~26 MB |
| STEP entities | 295,445 |
| Physical elements (`IfcElement`) | 17,553, all typed (65 type objects) |
| Georeferencing | `IfcMapConversion` to `EPSG:6594`, NAD83(2011) / Virginia South (m), NAVD88 heights |
| Alignments | 2 (`ALN-CBBT-US13-NB` ≈ 31.0 km, `ALN-CBBT-US13-SB` ≈ 31.2 km) |
| Referents | 20 `IfcReferent` station markers (STA 0+00 to 1011+77) |
| Classification | CCI Construction, US Bridge Data Dictionary, RIH Materials |

## Spatial structure

```
IfcProject
└─ IfcSite  Chesapeake Bay Mouth Corridor (Virginia Beach to Wise Point)
   ├─ IfcBridge          Trestles A–F, North Channel & Fisherman Inlet bridges   (8 IfcBridgePart)
   ├─ IfcRoad            US-13 approaches, toll plazas, island transitions, causeway   (4 IfcRoadPart)
   ├─ IfcMarineFacility  Portal Islands #1–#4, Sea Gull Pier, dolphins, nav beacons   (8 IfcMarinePart)
   └─ IfcFacility        Subaqueous tunnels: Thimble Shoal & Chesapeake Channel   (4 IfcFacilityPartCommon)
```

IFC 4.3 has no `IfcTunnel` entity, so the tunnels are modelled as a generic `IfcFacility`.

## Contents

| Entity | Count | | Entity | Count |
| :--- | ---: | --- | :--- | ---: |
| IfcCovering (mostly road striping) | 3,224 | | IfcWall | 166 |
| IfcPile | 2,419 | | IfcBuiltElement | 89 |
| IfcBeam | 2,221 | | IfcVehicle (67 road, 4 vessels) | 71 |
| IfcKerb | 1,958 | | IfcFooting | 60 |
| IfcRailing | 1,844 | | IfcMooringDevice | 30 |
| IfcColumn | 1,195 | | IfcNavigationElement | 17 |
| IfcPavement | 1,057 | | IfcAlarm / IfcSignal / IfcFireSuppressionTerminal / IfcAudioVisualAppliance | 14 each |
| IfcSlab | 865 | | IfcSensor | 8 |
| IfcLightFixture | 783 | | IfcUnitaryEquipment | 4 |
| IfcCourse | 571 | | IfcAlignmentSegment | 600 |
| IfcDiscreteAccessory | 484 | | | |
| IfcSign | 223 | | | |
| IfcMember | 208 | | | |

Modelled features include:

- the trestles with cylinder-pile bents and emergency pull-off bays (4th pile, extended cap, outer girder);
- the North Channel 1964 through-truss and 1999 plate-girder bridges on steel H-pile groups;
- the Fisherman Inlet bridge and its timber "ship funnel" fenders;
- the four portal islands, each with one ventilation building, wavewalls, open-cut ramps, ITS equipment and tow trucks;
- the immersed-tube tunnel segments with jet fans, standpipes and SOS boxes;
- both toll plazas.

## Properties and quantities

| Set | Attached to | Notes |
| :--- | :--- | :--- |
| `Pset_Condition`, `Pset_Warranty`, `Pset_ManufacturerTypeInformation` | 49 types → 17,414 occurrences | Placeholder values; identical condition grade and date on every type. |
| `usBrPset_BridgeIdentification`, `usBrPset_DesignLoading` | 17,553 / 17,414 occurrences | Placeholder values, assigned without regard to element or location (see known issues). |
| `Qto_*BaseQuantities` | most structural classes | Where present, volumes match the geometry (within 1% in a 1,500-element sample). |
| `Qto_BodyGeometryValidation` | 1,542 occurrences | GrossSurfaceArea / NetVolume for fixtures, equipment and vehicles. |

## Known issues

These issues are deliberately documented. They make the file a useful test case for audits that check whether data is *plausible*, not just whether it is *present*.

1. **Placeholder agency data.** `usBrPset_*` values cycle through five value sets regardless of where an element is. For example, tunnel elements are tagged with trestle maintenance sectors. Several DesignLoading fields hold unrelated text (e.g. `WindLoading` contains scour classes).
2. **Contradictory condition ratings.** `Pset_Condition` gives NBIS Grade 8 everywhere, while `TypeServiceOn` on the same elements gives Grade 7, 8 or 9.
3. **Partial quantities.** GrossSurfaceArea is on about 9% of elements and NetVolume on about 55%. Kerbs, railings/parapets, coverings and signs carry no volume.
4. **Type mismatches.**
   - Bridge parapets are typed as W-beam guardrail; the concrete parapet type has no occurrences.
   - The "2026 TBM tunnel lining" type is assigned to four 1964 immersed-tube elements.
   - Jet fans, floodgates, dolphins and toll booths are `IfcBuiltElement` instead of more specific classes (e.g. `IfcFan.JETFAN`).
5. **Not modelled.**
   - The 2026 Parallel Thimble Shoal Tunnel.
   - Individual armour units (hexapods). Each island's armour is a single element derived from an OSM shoreline polygon.
6. **Stationing is name-only.** All elements use `IfcLocalPlacement`; there is no `IfcLinearPlacement`. Referents carry a custom `CBBT_StationReferent` set rather than `Pset_Stationing`.
7. **Non-unique span numbers** in some pull-off bay element names (e.g. "Span 1" occurs several times on Trestle C NB).
8. **Asset count includes non-assets.** These include about 3,159 road-marking stripes, 67 passenger/freight vehicles and 4 vessels.

## Using the file

It opens in IFC 4.3-capable viewers (e.g. Bonsai/BlenderBIM, IfcOpenShell-based tools). It is large, so expect several seconds to load. Viewers without IFC 4.3 support will not show the facilities or alignments correctly.

## Sources and attribution

- **Geometry and layout:** approximated from public information and OpenStreetMap data (© OpenStreetMap contributors, ODbL).
- **Names:** facility and element names refer to the real CBBT for realism only.
- **Images:** aerial and Street View images used during development are not part of this repository.
