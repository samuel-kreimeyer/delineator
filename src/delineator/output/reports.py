"""Report generation for watershed analysis."""

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import geopandas as gpd

from delineator.config import ReportFormat
from delineator.exceptions import OutputError
from delineator.utils.logging import get_logger


class ReportGenerator:
    """Generates reports from watershed analysis results."""

    def __init__(self):
        """Initialize the report generator."""
        self.logger = get_logger(__name__)

    def generate(
        self,
        watersheds_gdf: gpd.GeoDataFrame,
        output_path: Path,
        report_format: ReportFormat,
    ) -> None:
        """Generate a report in the specified format.

        Args:
            watersheds_gdf: GeoDataFrame with watershed data.
            output_path: Path for output report file.
            report_format: Report format (json, csv, text).
        """
        self.logger.info(f"Generating {report_format.value} report")

        if report_format == ReportFormat.JSON:
            self._generate_json(watersheds_gdf, output_path)
        elif report_format == ReportFormat.CSV:
            self._generate_csv(watersheds_gdf, output_path)
        elif report_format == ReportFormat.TEXT:
            self._generate_text(watersheds_gdf, output_path)
        else:
            raise OutputError(f"Unsupported report format: {report_format}")

    def _generate_json(
        self, watersheds_gdf: gpd.GeoDataFrame, output_path: Path
    ) -> None:
        """Generate JSON report."""
        report_data = self._build_report_data(watersheds_gdf)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(report_data, f, indent=2, default=str)

    def _generate_csv(
        self, watersheds_gdf: gpd.GeoDataFrame, output_path: Path
    ) -> None:
        """Generate CSV report."""
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Select columns for CSV (exclude geometry) - imperial units
        columns = [
            "watershed_id",
            "point_name",
            "area_sqft",
            "area_acres",
            "area_sqmi",
            "elev_min",
            "elev_max",
            "elev_mean",
            "lfp_length",
            "lfp_slope",
            "lfp_elev_start",
            "lfp_elev_end",
            "lfp_elev_drop",
        ]
        available_columns = [c for c in columns if c in watersheds_gdf.columns]

        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(available_columns)

            for _, row in watersheds_gdf.iterrows():
                writer.writerow([row[c] for c in available_columns])

    def _generate_text(
        self, watersheds_gdf: gpd.GeoDataFrame, output_path: Path
    ) -> None:
        """Generate plain text report."""
        report_data = self._build_report_data(watersheds_gdf)

        output_path.parent.mkdir(parents=True, exist_ok=True)

        lines = [
            "=" * 60,
            "WATERSHED DELINEATION REPORT",
            "=" * 60,
            "",
            f"Generated: {report_data['metadata']['generated_at']}",
            f"Total watersheds: {report_data['metadata']['total_watersheds']}",
            "",
            "-" * 60,
            "SUMMARY STATISTICS",
            "-" * 60,
            "",
        ]

        summary = report_data["summary"]
        lines.extend([
            f"Total area: {summary['total_area_sqmi']:.4f} mi²",
            f"           {summary['total_area_acres']:.2f} acres",
            f"           {summary['total_area_sqft']:.0f} ft²",
            "",
            f"Elevation range: {summary['min_elevation_ft']:.2f} - {summary['max_elevation_ft']:.2f} ft",
            "",
            "-" * 60,
            "INDIVIDUAL WATERSHEDS",
            "-" * 60,
            "",
        ])

        for ws in report_data["watersheds"]:
            lines.extend([
                f"Watershed {ws['watershed_id']}: {ws['point_name']}",
                f"  Area: {ws['area_acres']:.2f} acres ({ws['area_sqmi']:.4f} mi²)",
                f"  Elevation: min={ws['elev_min_ft']:.2f} ft, max={ws['elev_max_ft']:.2f} ft, mean={ws['elev_mean_ft']:.2f} ft",
            ])

            # Hydraulic distance (longest flow path) information
            if ws.get("hydraulic_distance_ft") is not None:
                lines.append(f"  Hydraulic distance (longest flow path): {ws['hydraulic_distance_ft']:.2f} ft")
                if ws.get("lfp_slope") is not None:
                    lines.append(f"    Average slope: {ws['lfp_slope']:.6f} ft/ft ({ws['lfp_slope']*100:.4f}%)")
                if ws.get("lfp_elev_drop_ft") is not None:
                    lines.append(f"    Elevation drop: {ws['lfp_elev_drop_ft']:.2f} ft")
                if ws.get("lfp_elev_start_ft") is not None and ws.get("lfp_elev_end_ft") is not None:
                    lines.append(f"    Start elevation: {ws['lfp_elev_start_ft']:.2f} ft, End elevation: {ws['lfp_elev_end_ft']:.2f} ft")
            else:
                lines.append("  Hydraulic distance (longest flow path): N/A")

            lines.append("")

        lines.extend([
            "=" * 60,
            "END OF REPORT",
            "=" * 60,
        ])

        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

    def _build_report_data(self, watersheds_gdf: gpd.GeoDataFrame) -> dict[str, Any]:
        """Build report data structure with imperial units."""
        watersheds = []

        for _, row in watersheds_gdf.iterrows():
            ws_data = {
                "watershed_id": int(row["watershed_id"]),
                "point_name": row.get("point_name", f"Point_{row['watershed_id']}"),
                "area_sqft": float(row.get("area_sqft", 0)),
                "area_acres": float(row.get("area_acres", 0)),
                "area_sqmi": float(row.get("area_sqmi", 0)),
                "elev_min_ft": float(row.get("elev_min", 0)) if row.get("elev_min") is not None else None,
                "elev_max_ft": float(row.get("elev_max", 0)) if row.get("elev_max") is not None else None,
                "elev_mean_ft": float(row.get("elev_mean", 0)) if row.get("elev_mean") is not None else None,
                # Longest flow path (hydraulic distance) in feet (if available)
                "hydraulic_distance_ft": float(row.get("lfp_length")) if row.get("lfp_length") is not None else None,
                "lfp_slope": float(row.get("lfp_slope")) if row.get("lfp_slope") is not None else None,
                "lfp_elev_start_ft": float(row.get("lfp_elev_start")) if row.get("lfp_elev_start") is not None else None,
                "lfp_elev_end_ft": float(row.get("lfp_elev_end")) if row.get("lfp_elev_end") is not None else None,
                "lfp_elev_drop_ft": float(row.get("lfp_elev_drop")) if row.get("lfp_elev_drop") is not None else None,
            }
            watersheds.append(ws_data)

        # Compute summary
        total_area_sqft = sum(ws["area_sqft"] for ws in watersheds)
        total_area_acres = sum(ws["area_acres"] for ws in watersheds)
        total_area_sqmi = sum(ws["area_sqmi"] for ws in watersheds)

        elevations = [
            e for ws in watersheds
            for e in [ws["elev_min_ft"], ws["elev_max_ft"]]
            if e is not None
        ]

        summary = {
            "total_area_sqft": total_area_sqft,
            "total_area_acres": total_area_acres,
            "total_area_sqmi": total_area_sqmi,
            "min_elevation_ft": min(elevations) if elevations else None,
            "max_elevation_ft": max(elevations) if elevations else None,
        }

        return {
            "metadata": {
                "generated_at": datetime.now().isoformat(),
                "total_watersheds": len(watersheds),
            },
            "summary": summary,
            "watersheds": watersheds,
        }
