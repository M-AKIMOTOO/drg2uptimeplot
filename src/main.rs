use astro::{coords, ecliptic, sun, time};
use chrono::{DateTime, Datelike, Duration, NaiveDate, NaiveDateTime, TimeZone, Timelike, Utc};
use clap::Parser;
use eframe::egui;
use egui::{Color32, Stroke};
use egui_plot::{Corner, GridMark, Legend, Line, LineStyle, Plot, PlotBounds, PlotPoints, Points};
use nav_types::{ECEF, WGS84};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

// --- Command Line Arguments Definition ---
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Path to the DRG schedule file
    drg_file: String,

    /// Name of the station to use for calculations
    #[arg(long, default_value = "YAMAGU32")]
    station: String,

    /// Print the AZ/EL samples used by the GUI and exit
    #[arg(long)]
    terminal: bool,

    /// Write full sky tracks with elevation >= 0 degrees
    #[arg(long, value_name = "FILE", num_args = 0..=1, default_missing_value = "DRG_azel.txt")]
    output: Option<String>,
}

// --- Data Structures ---
#[derive(Debug, Clone)]
struct Source {
    name1: String,
    name2: String,
    ra_rad: f64,
    dec_rad: f64,
}

#[derive(Debug, Clone)]
struct Station {
    name: String,
    pos: [f64; 3],
}

#[derive(Debug, Clone)]
struct Observation {
    source_name: String,
    start_time: NaiveDateTime,
    duration_sec: i64,
}

#[derive(Debug)]
struct DrgData {
    sources: Vec<Source>,
    schedule: Vec<Observation>,
}

#[derive(PartialEq)]
enum AppTab {
    UptimePlot,
    UptimePlot2,
    PolarPlot1,
    PolarPlot2,
}

// --- Plotting App ---
struct DrgPlotApp {
    sky_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
    drg_ut_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
    plot_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
    color_map: HashMap<String, Color32>,
    x_axis_bounds: [f64; 2],
    t0: NaiveDateTime,
    selected_tab: AppTab,
    reset_plot_bounds: bool,
}

impl DrgPlotApp {
    fn new(
        station: Station,
        drg_data: DrgData,
        x_axis_bounds: [f64; 2],
        t0: NaiveDateTime,
        t_end: NaiveDateTime,
    ) -> Self {
        let sky_segments =
            calculate_ut_sky_segments(&station, &drg_data.sources, &drg_data.schedule, t0, t_end);
        let drg_ut_segments =
            calculate_observation_segments_ut(&station, &drg_data.sources, &drg_data.schedule);
        let mut plot_segments =
            calculate_observation_segments(&station, &drg_data.sources, &drg_data.schedule, t0);
        let sun_segment = calculate_sun_segments(&station, t0, t_end, false);
        plot_segments.push(sun_segment);

        let mut color_map = HashMap::new();
        let palette = [
            Color32::from_rgb(100, 143, 255), // Blue
            Color32::from_rgb(120, 255, 120), // Green
            Color32::from_rgb(255, 100, 100), // Red
            Color32::from_rgb(255, 180, 80),  // Orange
            Color32::from_rgb(240, 120, 240), // Magenta
            Color32::from_rgb(130, 255, 255), // Cyan
        ];

        let mut unique_sources = drg_data
            .sources
            .iter()
            .map(|s| s.name2.clone())
            .collect::<Vec<_>>();
        unique_sources.sort();
        unique_sources.dedup();

        for (i, name) in unique_sources.iter().enumerate() {
            color_map.insert(name.clone(), palette[i % palette.len()]);
        }
        color_map.insert("Sun".to_string(), Color32::from_rgb(255, 255, 0)); // Yellow for Sun

        Self {
            sky_segments,
            drg_ut_segments,
            plot_segments,
            color_map,
            x_axis_bounds,
            t0,
            selected_tab: AppTab::UptimePlot,
            reset_plot_bounds: false,
        }
    }
}

impl eframe::App for DrgPlotApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.selected_tab, AppTab::UptimePlot, "UptimePlot1");
                ui.selectable_value(&mut self.selected_tab, AppTab::UptimePlot2, "UptimePlot2");
                ui.selectable_value(&mut self.selected_tab, AppTab::PolarPlot1, "PolarPlot1");
                ui.selectable_value(&mut self.selected_tab, AppTab::PolarPlot2, "PolarPlot2");

                ui.separator();

                if ui.button("Reset Zoom").clicked() {
                    self.reset_plot_bounds = true;
                }
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| match self.selected_tab {
            AppTab::UptimePlot => self.ui_uptime_plot_tab(ui),
            AppTab::UptimePlot2 => self.ui_uptime_plot2_tab(ui),
            AppTab::PolarPlot1 => self.ui_polar_plot1_tab(ui),
            AppTab::PolarPlot2 => self.ui_polar_plot2_tab(ui),
        });
        self.reset_plot_bounds = false;
    }
}

impl DrgPlotApp {
    fn ui_uptime_plot_tab(&mut self, ui: &mut egui::Ui) {
        let t0 = self.t0;
        let pointer_time_formatter = move |x: f64| -> String {
            let time = t0 + Duration::seconds((x * 3600.0).round() as i64);
            time.format("%m-%d %H:%M").to_string()
        };

        let az_pointer_formatter = |x: f64, y: f64| {
            format!(
                "Time: {:02}:{:02}\nAz: {:.1}°",
                x as u32,
                (x.fract() * 60.0) as u32,
                y
            )
        };
        let el_pointer_formatter = |x: f64, y: f64| {
            format!(
                "Time: {:02}:{:02}\nEl: {:.1}°",
                x as u32,
                (x.fract() * 60.0) as u32,
                y
            )
        };

        let plot_az = Plot::new("az_plot")
            .width(ui.available_width())
            .height(ui.available_height() / 2.0)
            .y_axis_label("Azimuth (deg)")
            .y_axis_min_width(70.0) // Changed from 0.0 to 70.0 for uniform width
            .allow_drag(true)
            .allow_zoom(true)
            .allow_scroll(true)
            .include_y(0.0)
            .include_y(360.0)
            .y_grid_spacer(|_input| {
                (0..=12)
                    .map(|v| GridMark {
                        value: (v * 30) as f64,
                        step_size: 30.0,
                    })
                    .collect()
            }) // Re-added
            .show_x(true) // Re-added, Temporarily set to true for testing
            .x_axis_label("") // Re-added
            .x_axis_formatter(|_, _| "".to_string()) // Re-added
            .y_axis_formatter(|m, _| format!("{:>3}", m.value as i32))
            .show_y(true) // Re-added, changed to 3-digit padded
            .coordinates_formatter(
                Corner::LeftTop,
                egui_plot::CoordinatesFormatter::new(move |p, _| az_pointer_formatter(p.x, p.y)),
            )
            .legend(Legend::default());

        let plot_el = Plot::new("el_plot")
            .width(ui.available_width())
            .height(ui.available_height() / 2.0)
            .y_axis_label("Elevation (deg)")
            .y_axis_min_width(70.0) // Changed from 67.0 to 70.0 for uniform width
            .allow_drag(true)
            .allow_zoom(true)
            .allow_scroll(true)
            .include_y(0.0)
            .include_y(90.0)
            .x_axis_formatter(move |m, _| pointer_time_formatter(m.value)) // Re-added
            .y_axis_formatter(|m, _| format!("{:>2}", m.value as i32))
            .show_y(true) // Re-added, changed to 3-digit padded
            .y_grid_spacer(|_input| {
                (0..=9)
                    .map(|v| GridMark {
                        value: (v * 10) as f64,
                        step_size: 10.0,
                    })
                    .collect()
            }) // Re-added
            .coordinates_formatter(
                Corner::LeftTop,
                egui_plot::CoordinatesFormatter::new(move |p, _| el_pointer_formatter(p.x, p.y)),
            )
            .legend(Legend::default());

        plot_az.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max(
                    [self.x_axis_bounds[0], 0.0],
                    [self.x_axis_bounds[1], 360.0],
                ));
            }
            for (name, az_points, _) in &self.plot_segments {
                if let Some(color) = self.color_map.get(name) {
                    plot_ui.line(
                        Line::new(name.clone(), PlotPoints::from(az_points.clone())).color(*color),
                    );
                }
            }
        });

        ui.add_space(-10.0);

        plot_el.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max(
                    [self.x_axis_bounds[0], 0.0],
                    [self.x_axis_bounds[1], 90.0],
                ));
            }
            for (name, _, el_points) in &self.plot_segments {
                if let Some(color) = self.color_map.get(name) {
                    plot_ui.line(
                        Line::new(name.clone(), PlotPoints::from(el_points.clone())).color(*color),
                    );
                }
            }
        });
    }

    fn ui_uptime_plot2_tab(&mut self, ui: &mut egui::Ui) {
        let pointer_time_formatter = |x: f64| -> String {
            let total_minutes = (x * 60.0).round() as i64;
            format!("{:02}:{:02}", total_minutes / 60, total_minutes % 60)
        };

        let az_pointer_formatter = |x: f64, y: f64| {
            format!(
                "Time: {:02}:{:02}\nAz: {:.1}°",
                x as u32,
                (x.fract() * 60.0) as u32,
                y
            )
        };
        let el_pointer_formatter = |x: f64, y: f64| {
            format!(
                "Time: {:02}:{:02}\nEl: {:.1}°",
                x as u32,
                (x.fract() * 60.0) as u32,
                y
            )
        };

        let plot_az = Plot::new("az_plot_ut24")
            .width(ui.available_width())
            .default_x_bounds(0.0, 24.0)
            .height(ui.available_height() / 2.0)
            .y_axis_label("Azimuth (deg)")
            .y_axis_min_width(70.0)
            .allow_drag(true)
            .allow_zoom(true)
            .allow_scroll(true)
            .include_y(0.0)
            .include_y(360.0)
            .y_grid_spacer(|_input| {
                (0..=12)
                    .map(|v| GridMark {
                        value: (v * 30) as f64,
                        step_size: 30.0,
                    })
                    .collect()
            })
            .show_x(true)
            .x_axis_label("")
            .x_axis_formatter(|_, _| "".to_string())
            .y_axis_formatter(|m, _| format!("{:>3}", m.value as i32))
            .show_y(true)
            .coordinates_formatter(
                Corner::LeftTop,
                egui_plot::CoordinatesFormatter::new(move |p, _| az_pointer_formatter(p.x, p.y)),
            )
            .legend(Legend::default());

        let plot_el = Plot::new("el_plot_ut24")
            .width(ui.available_width())
            .default_x_bounds(0.0, 24.0)
            .height(ui.available_height() / 2.0)
            .y_axis_label("Elevation (deg)")
            .y_axis_min_width(70.0)
            .allow_drag(true)
            .allow_zoom(true)
            .allow_scroll(true)
            .default_y_bounds(0.0, 90.0)
            .include_y(0.0)
            .include_y(90.0)
            .x_axis_label("UTC time")
            .x_axis_formatter(move |m, _| pointer_time_formatter(m.value))
            .y_axis_formatter(|m, _| format!("{:>2}", m.value as i32))
            .show_y(true)
            .y_grid_spacer(|_input| {
                (0..=9)
                    .map(|v| GridMark {
                        value: (v * 10) as f64,
                        step_size: 10.0,
                    })
                    .collect()
            })
            .coordinates_formatter(
                Corner::LeftTop,
                egui_plot::CoordinatesFormatter::new(move |p, _| el_pointer_formatter(p.x, p.y)),
            )
            .legend(Legend::default());

        plot_az.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([0.0, 0.0], [24.0, 360.0]));
            }

            for (name, az_points, _) in &self.sky_segments {
                if let Some(color) = self.color_map.get(name) {
                    plot_ui.line(
                        Line::new(
                            format!("{} (sky)", name),
                            PlotPoints::from(az_points.clone()),
                        )
                        .color(*color)
                        .width(1.0)
                        .style(LineStyle::dotted_dense()),
                    );
                }
            }

            let mut seen_drg_targets = HashSet::new();
            for (name, az_points, _) in &self.drg_ut_segments {
                if name == "Sun" {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let label = if seen_drg_targets.insert(name.clone()) {
                        format!("{} (DRG)", name)
                    } else {
                        String::new()
                    };
                    plot_ui.line(
                        Line::new(label, PlotPoints::from(az_points.clone()))
                            .color(*color)
                            .width(3.0),
                    );
                }
            }
        });

        ui.add_space(-10.0);

        plot_el.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([0.0, 0.0], [24.0, 90.0]));
            }

            for (name, _, el_points) in &self.sky_segments {
                if let Some(color) = self.color_map.get(name) {
                    let above_horizon = el_points
                        .iter()
                        .copied()
                        .filter(|point| point[1] > 0.0)
                        .collect::<Vec<_>>();
                    plot_ui.line(
                        Line::new(
                            format!("{} (sky)", name),
                            PlotPoints::from(above_horizon),
                        )
                        .color(*color)
                        .width(1.0)
                        .style(LineStyle::dotted_dense()),
                    );
                }
            }

            let mut seen_drg_targets = HashSet::new();
            for (name, _, el_points) in &self.drg_ut_segments {
                if name == "Sun" {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let label = if seen_drg_targets.insert(name.clone()) {
                        format!("{} (DRG)", name)
                    } else {
                        String::new()
                    };
                    let above_horizon = el_points
                        .iter()
                        .copied()
                        .filter(|point| point[1] > 0.0)
                        .collect::<Vec<_>>();
                    plot_ui.line(
                        Line::new(label, PlotPoints::from(above_horizon))
                            .color(*color)
                            .width(3.0),
                    );
                }
            }
        });
    }

    fn ui_polar_plot1_tab(&mut self, ui: &mut egui::Ui) {
        let plot = Plot::new("polar_plot")
            .width(ui.available_width())
            .height(ui.available_height())
            .data_aspect(1.0)
            .view_aspect(1.0)
            .include_x(-1.0)
            .include_x(1.0)
            .include_y(-1.0)
            .include_y(1.0)
            .center_x_axis(true)
            .center_y_axis(true)
            .show_x(false)
            .show_y(false)
            .x_grid_spacer(|_input| vec![])
            .y_grid_spacer(|_input| vec![])
            .allow_drag(true)
            .allow_zoom(true)
            .allow_scroll(true) // Added for interactivity
            .legend(Legend::default());

        plot.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([-1.0, -1.0], [1.0, 1.0]));
            }
            draw_polar_grid(plot_ui, 0.0);

            for (name, az_points, el_points) in &self.plot_segments {
                let mut polar_points = Vec::new();
                for i in 0..az_points.len() {
                    let az = az_points[i][1];
                    let el = el_points[i][1];
                    if !el.is_nan() && el >= 0.0 {
                        let angle_rad = (90.0f64 - az).to_radians();
                        let radius = (90.0 - el) / 90.0;
                        polar_points.push([radius * angle_rad.cos(), radius * angle_rad.sin()]);
                    }
                }
                if !polar_points.is_empty() {
                    if let Some(color) = self.color_map.get(name) {
                        plot_ui.points(
                            Points::new(name.clone(), PlotPoints::from(polar_points)).color(*color),
                        );
                    }
                }
            }
        });
    }

    fn ui_polar_plot2_tab(&mut self, ui: &mut egui::Ui) {
        let plot = Plot::new("polar_plot2")
            .width(ui.available_width())
            .height(ui.available_height())
            .data_aspect(1.0)
            .view_aspect(1.0)
            .include_x(-1.0)
            .include_x(1.0)
            .include_y(-1.0)
            .include_y(1.0)
            .center_x_axis(true)
            .center_y_axis(true)
            .show_x(false)
            .show_y(false)
            .x_grid_spacer(|_input| vec![])
            .y_grid_spacer(|_input| vec![])
            .allow_drag(true)
            .allow_zoom(true)
            .allow_scroll(true)
            .legend(Legend::default());

        plot.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([-1.0, -1.0], [1.0, 1.0]));
            }
            draw_polar_grid(plot_ui, 0.0);

            for (name, az_points, el_points) in &self.sky_segments {
                if let Some(color) = self.color_map.get(name) {
                    let points = azel_to_polar_points(az_points, el_points);
                    plot_ui.line(
                        Line::new(format!("{} (sky)", name), PlotPoints::from(points))
                            .color(*color)
                            .width(1.0)
                            .style(LineStyle::dotted_dense()),
                    );
                }
            }

            let mut seen_drg_targets = HashSet::new();
            for (name, az_points, el_points) in &self.drg_ut_segments {
                if name == "Sun" {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let label = if seen_drg_targets.insert(name.clone()) {
                        format!("{} (DRG)", name)
                    } else {
                        String::new()
                    };
                    let points = azel_to_polar_points(az_points, el_points);
                    plot_ui.line(
                        Line::new(label, PlotPoints::from(points))
                            .color(*color)
                            .width(3.0),
                    );
                }
            }
        });
    }
}

fn draw_polar_grid(plot_ui: &mut egui_plot::PlotUi<'_>, minimum_elevation: f64) {
    let ring_count = ((90.0 - minimum_elevation) / 15.0).round() as usize;
    for ring in 0..=ring_count {
        let elevation = minimum_elevation + ring as f64 * 15.0;
        let radius = (90.0 - elevation) / 90.0;
        let circle_points: PlotPoints = (0..=100)
            .map(|i| {
                let angle = i as f64 * 2.0 * std::f64::consts::PI / 100.0;
                [radius * angle.cos(), radius * angle.sin()]
            })
            .collect();
        plot_ui.line(Line::new("", circle_points).stroke(Stroke::new(1.0, Color32::DARK_GRAY)));
        if elevation != 90.0 {
            let label_angle = 72.0f64.to_radians();
            plot_ui.text(
                egui_plot::Text::new(
                    "",
                    egui_plot::PlotPoint::new(
                        radius * label_angle.cos(),
                        radius * label_angle.sin(),
                    ),
                    format!("{:.0}°", elevation),
                )
                .color(Color32::DARK_GRAY),
            );
        }
    }

    let maximum_radius = (90.0 - minimum_elevation) / 90.0;
    for azimuth in [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0] {
        let angle = (90.0f64 - azimuth).to_radians();
        let line_points = vec![
            [0.0, 0.0],
            [maximum_radius * angle.cos(), maximum_radius * angle.sin()],
        ];
        plot_ui.line(
            Line::new("", PlotPoints::from(line_points))
                .stroke(Stroke::new(1.0, Color32::DARK_GRAY)),
        );
        plot_ui.text(
            egui_plot::Text::new(
                "",
                egui_plot::PlotPoint::new(
                    (maximum_radius + 0.1) * angle.cos(),
                    (maximum_radius + 0.1) * angle.sin(),
                ),
                format!("{:.0}°", azimuth),
            )
            .color(Color32::DARK_GRAY),
        );
    }
}

fn azel_to_polar_points(az_points: &[[f64; 2]], el_points: &[[f64; 2]]) -> Vec<[f64; 2]> {
    az_points
        .iter()
        .zip(el_points)
        .map(|(az_point, el_point)| {
            let azimuth = az_point[1];
            let elevation = el_point[1];
            if !azimuth.is_finite() || !elevation.is_finite() || elevation <= 0.0 {
                return [f64::NAN, f64::NAN];
            }
            let angle = (90.0 - azimuth).to_radians();
            let radius = (90.0 - elevation) / 90.0;
            [radius * angle.cos(), radius * angle.sin()]
        })
        .collect()
}

pub fn radec2azalt(
    ant_position: [f32; 3],
    time: DateTime<Utc>,
    obs_ra: f32,
    obs_dec: f32,
) -> (f32, f32, f32) {
    let obs_year = time.year() as i16;
    let obs_month = time.month() as u8;
    let obs_day = time.day() as u8;
    let obs_hour = time.hour() as u8;
    let obs_minute = time.minute() as u8;
    let obs_second = time.second() as f64; // + (time.nanosecond() as f64 / 1_000_000_000.0);

    let decimal_day_calc = obs_day as f64
        + obs_hour as f64 / 24.0
        + obs_minute as f64 / 60.0 / 24.0
        + obs_second as f64 / 24.0 / 60.0 / 60.0;

    let date = time::Date {
        year: obs_year,
        month: obs_month,
        decimal_day: decimal_day_calc,
        cal_type: time::CalType::Gregorian,
    };

    let ecef_position = ECEF::new(
        ant_position[0] as f64,
        ant_position[1] as f64,
        ant_position[2] as f64,
    );
    let wgs84_position: WGS84<f64> = ecef_position.into();
    let longitude_radian = wgs84_position.longitude_radians();
    let latitude_radian = wgs84_position.latitude_radians();
    let height_meter = wgs84_position.altitude();

    let julian_day = time::julian_day(&date);
    let mean_sidereal = time::mn_sidr(julian_day);
    let hour_angle =
        coords::hr_angl_frm_observer_long(mean_sidereal, -longitude_radian, obs_ra as f64);

    let source_az =
        coords::az_frm_eq(hour_angle, obs_dec as f64, latitude_radian).to_degrees() as f32 + 180.0;
    let source_el =
        coords::alt_frm_eq(hour_angle, obs_dec as f64, latitude_radian).to_degrees() as f32;

    (source_az, source_el, height_meter as f32)
}

// --- Calculation Logic ---

fn calculate_sun_segments(
    station: &Station,
    t0: NaiveDateTime,
    t_end: NaiveDateTime,
    include_below_horizon: bool,
) -> (String, Vec<[f64; 2]>, Vec<[f64; 2]>) {
    let mut az_segment = Vec::new();
    let mut el_segment = Vec::new();
    let ant_pos = station.pos;

    let mut current_time = t0;
    while current_time <= t_end {
        let duration_since_t0 = current_time.signed_duration_since(t0);
        let hour_float = duration_since_t0.num_seconds() as f64 / 3600.0;

        let datetime_utc = Utc.from_utc_datetime(&current_time);

        let obs_year = datetime_utc.year() as i16;
        let obs_month = datetime_utc.month() as u8;
        let obs_day = datetime_utc.day() as u8;
        let obs_hour = datetime_utc.hour() as u8;
        let obs_minute = datetime_utc.minute() as u8;
        let obs_second = datetime_utc.second() as f64;
        let decimal_day_calc = obs_day as f64
            + obs_hour as f64 / 24.0
            + obs_minute as f64 / 60.0 / 24.0
            + obs_second as f64 / 24.0 / 60.0 / 60.0;

        let date = time::Date {
            year: obs_year,
            month: obs_month,
            decimal_day: decimal_day_calc,
            cal_type: time::CalType::Gregorian,
        };
        let jd = time::julian_day(&date);

        // 1. Get Sun's ecliptic coordinates
        let (ecl_point, _) = sun::geocent_ecl_pos(jd);
        let ecl_lon = ecl_point.long;
        let ecl_lat = ecl_point.lat;

        // 2. Convert ecliptic to equatorial (RA/Dec)
        let ra_rad = coords::asc_frm_ecl(ecl_lon, ecl_lat, ecliptic::mn_oblq_IAU(jd));
        let dec_rad = coords::dec_frm_ecl(ecl_lon, ecl_lat, ecliptic::mn_oblq_IAU(jd));

        // 3. Convert equatorial (RA/Dec) to Az/El
        let mean_sidereal = time::mn_sidr(jd);
        let geocentric_coord = ECEF::new(ant_pos[0] as f64, ant_pos[1] as f64, ant_pos[2] as f64);
        let geodetic_coord: WGS84<f64> = geocentric_coord.into();
        let longitude_radian = geodetic_coord.longitude_radians();
        let latitude_radian = geodetic_coord.latitude_radians();
        let hour_angle =
            coords::hr_angl_frm_observer_long(mean_sidereal, -longitude_radian, ra_rad);

        let az =
            coords::az_frm_eq(hour_angle, dec_rad, latitude_radian).to_degrees() as f32 + 180.0;
        let el = coords::alt_frm_eq(hour_angle, dec_rad, latitude_radian).to_degrees() as f32;

        if el >= 0.0 || include_below_horizon {
            az_segment.push([hour_float, az as f64]);
            el_segment.push([hour_float, el as f64]);
        } else {
            az_segment.push([hour_float, az as f64]);
            el_segment.push([hour_float, f64::NAN]);
        }
        current_time += Duration::minutes(5); // Calculate every 5 minutes
    }

    ("Sun".to_string(), az_segment, el_segment)
}

fn calculate_sky_segments(
    station: &Station,
    sources: &[Source],
    schedule: &[Observation],
    t0: NaiveDateTime,
    t_end: NaiveDateTime,
) -> Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)> {
    let mut sky_segments = Vec::new();
    let mut seen_sources = HashSet::new();
    let scheduled_sources = schedule
        .iter()
        .map(|observation| observation.source_name.as_str())
        .collect::<HashSet<_>>();

    for source in sources {
        if !scheduled_sources.contains(source.name1.as_str())
            && !scheduled_sources.contains(source.name2.as_str())
        {
            continue;
        }
        if !seen_sources.insert(source.name2.clone()) {
            continue;
        }

        let mut az_segment = Vec::new();
        let mut el_segment = Vec::new();
        let mut current_time = t0;

        while current_time <= t_end {
            let duration_since_t0 = current_time.signed_duration_since(t0);
            let hour_float = duration_since_t0.num_seconds() as f64 / 3600.0;
            let datetime_utc = Utc.from_utc_datetime(&current_time);
            let (az, el, _) = radec2azalt(
                [
                    station.pos[0] as f32,
                    station.pos[1] as f32,
                    station.pos[2] as f32,
                ],
                datetime_utc,
                source.ra_rad as f32,
                source.dec_rad as f32,
            );

            az_segment.push([hour_float, az as f64]);
            el_segment.push([hour_float, el as f64]);
            current_time += Duration::minutes(1);
        }

        sky_segments.push((source.name2.clone(), az_segment, el_segment));
    }

    sky_segments.push(calculate_sun_segments(station, t0, t_end, true));
    sky_segments
}

fn utc_hour(time: NaiveDateTime) -> f64 {
    time.hour() as f64 + time.minute() as f64 / 60.0 + time.second() as f64 / 3600.0
}

fn convert_segments_to_ut_axis(
    base_time: NaiveDateTime,
    segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
) -> Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)> {
    segments
        .into_iter()
        .map(|(name, az_points, el_points)| {
            let mut ut_az_points = Vec::new();
            let mut ut_el_points = Vec::new();
            let mut previous_date = None;

            for (az_point, el_point) in az_points.iter().zip(el_points.iter()) {
                let sample_time =
                    base_time + Duration::seconds((az_point[0] * 3600.0).round() as i64);
                if previous_date.is_some() && previous_date != Some(sample_time.date()) {
                    ut_az_points.push([0.0, f64::NAN]);
                    ut_el_points.push([0.0, f64::NAN]);
                }
                ut_az_points.push([utc_hour(sample_time), az_point[1]]);
                ut_el_points.push([utc_hour(sample_time), el_point[1]]);
                previous_date = Some(sample_time.date());
            }

            (name, ut_az_points, ut_el_points)
        })
        .collect()
}

fn calculate_ut_sky_segments(
    station: &Station,
    sources: &[Source],
    schedule: &[Observation],
    t0: NaiveDateTime,
    t_end: NaiveDateTime,
) -> Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)> {
    let first_day_start = t0.date().and_hms_opt(0, 0, 0).unwrap();
    let last_day_end = t_end.date().and_hms_opt(23, 59, 0).unwrap();
    let segments =
        calculate_sky_segments(station, sources, schedule, first_day_start, last_day_end);
    convert_segments_to_ut_axis(first_day_start, segments)
}

fn calculate_observation_segments_ut(
    station: &Station,
    sources: &[Source],
    schedule: &[Observation],
) -> Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)> {
    let t0 = schedule
        .iter()
        .map(|observation| observation.start_time)
        .min()
        .unwrap();
    convert_segments_to_ut_axis(
        t0,
        calculate_observation_segments(station, sources, schedule, t0),
    )
}

fn calculate_observation_segments(
    station: &Station,
    sources: &[Source],
    schedule: &[Observation],
    t0: NaiveDateTime,
) -> Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)> {
    let mut new_plot_data = Vec::new();
    let ant_pos = station.pos;

    for obs in schedule {
        if let Some(source) = sources
            .iter()
            .find(|s| s.name1 == obs.source_name || s.name2 == obs.source_name)
        {
            let mut az_segment = Vec::new();
            let mut el_segment = Vec::new();

            let start_time = obs.start_time;
            let end_time = start_time + Duration::seconds(obs.duration_sec);

            let mut current_time = start_time;
            while current_time <= end_time {
                let duration_since_t0 = current_time.signed_duration_since(t0);
                let hour_float = duration_since_t0.num_seconds() as f64 / 3600.0;

                let datetime_utc = Utc.from_utc_datetime(&current_time);
                let (az, el, _) = radec2azalt(
                    [ant_pos[0] as f32, ant_pos[1] as f32, ant_pos[2] as f32],
                    datetime_utc,
                    source.ra_rad as f32,
                    source.dec_rad as f32,
                );

                if el >= 0.0 {
                    az_segment.push([hour_float, az as f64]);
                    el_segment.push([hour_float, el as f64]);
                } else {
                    az_segment.push([hour_float, az as f64]);
                    el_segment.push([hour_float, f64::NAN]);
                }
                current_time += Duration::minutes(1);
            }
            new_plot_data.push((source.name2.clone(), az_segment, el_segment));
        }
    }
    new_plot_data
}

fn write_full_track_file<P: AsRef<Path>>(
    path: P,
    station: &Station,
    sources: &[Source],
    schedule: &[Observation],
    t0: NaiveDateTime,
    t_end: NaiveDateTime,
) -> Result<usize, Box<dyn std::error::Error>> {
    let sky_segments = calculate_sky_segments(station, sources, schedule, t0, t_end);
    let mut file = File::create(path)?;
    writeln!(file, "# Full sky track (EL >= 0 deg)")?;
    writeln!(file, "# start: {}", t0.format("%Y-%m-%d %H:%M:%S"))?;
    writeln!(file, "# end:   {}", t_end.format("%Y-%m-%d %H:%M:%S"))?;
    writeln!(
        file,
        "{:<16} {:19} {:>8} {:>8}",
        "target", "time", "AZ", "EL"
    )?;

    let mut row_count = 0;
    for (name, az_points, el_points) in sky_segments {
        for (az_point, el_point) in az_points.iter().zip(el_points.iter()) {
            let el = el_point[1];
            if !el.is_finite() || el < 0.0 {
                continue;
            }

            let sample_time = t0 + Duration::seconds((az_point[0] * 3600.0).round() as i64);
            writeln!(
                file,
                "{:<16} {} {:>8.2} {:>8.2}",
                name,
                sample_time.format("%Y-%m-%d %H:%M:%S"),
                az_point[1],
                el
            )?;
            row_count += 1;
        }
        writeln!(file)?;
    }

    Ok(row_count)
}

fn print_terminal_observation_data(
    station: &Station,
    sources: &[Source],
    schedule: &[Observation],
) {
    println!(
        "{:<16} {:19} {:19} {:>6}",
        "target", "start", "end", "length"
    );

    for obs in schedule {
        let Some(source) = sources
            .iter()
            .find(|source| source.name1 == obs.source_name || source.name2 == obs.source_name)
        else {
            eprintln!("Warning: source '{}' was not found", obs.source_name);
            continue;
        };

        let end_time = obs.start_time + Duration::seconds(obs.duration_sec);
        println!(
            "{:<16} {} {} {:>6}",
            source.name2,
            obs.start_time.format("%Y-%m-%d %H:%M:%S"),
            end_time.format("%Y-%m-%d %H:%M:%S"),
            obs.duration_sec
        );
        println!("  {:19} {:>8} {:>8}", "time", "AZ", "EL");

        let mut sample_time = obs.start_time;
        while sample_time <= end_time {
            let datetime_utc = Utc.from_utc_datetime(&sample_time);
            let (az, el, _) = radec2azalt(
                [
                    station.pos[0] as f32,
                    station.pos[1] as f32,
                    station.pos[2] as f32,
                ],
                datetime_utc,
                source.ra_rad as f32,
                source.dec_rad as f32,
            );
            let el_text = if el >= 0.0 {
                format!("{:8.2}", el)
            } else {
                "      --".to_string()
            };
            println!(
                "  {} {:>8.2} {}",
                sample_time.format("%Y-%m-%d %H:%M:%S"),
                az,
                el_text
            );
            sample_time += Duration::minutes(1);
        }
        println!();
    }
}

// --- File Parsing Logic ---
fn parse_drg_file<P: AsRef<Path>>(path: P) -> Result<DrgData, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sources = Vec::new();
    let mut schedule = Vec::new();
    enum ParseSection {
        None,
        Sources,
        Sked,
    }
    let mut current_section = ParseSection::None;

    for line in reader.lines() {
        let line = line?;
        let trimmed_line = line.trim();
        if trimmed_line.starts_with('$') {
            current_section = match trimmed_line {
                "$SOURCES" => ParseSection::Sources,
                "$SKED" => ParseSection::Sked,
                _ => ParseSection::None,
            };
            continue;
        }
        if trimmed_line.is_empty() || trimmed_line.starts_with('*') {
            continue;
        }

        match current_section {
            ParseSection::Sources => {
                let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
                if parts.len() >= 9 && parts[8] == "2000.0" {
                    let mut name1 = parts[0].to_string();
                    let mut name2 = parts[1].to_string();
                    if name1 == "$" {
                        name1 = name2.clone();
                    }
                    if name2 == "$" {
                        name2 = name1.clone();
                    }
                    if name1 == "$" && name2 == "$" {
                        name1 = "NanashinoGonbei".to_string();
                        name2 = "NanashinoGonbei".to_string();
                    }
                    let ra_h: f64 = parts[2].parse()?;
                    let ra_m: f64 = parts[3].parse()?;
                    let ra_s: f64 = parts[4].parse()?;
                    let ra_hours = ra_h + ra_m / 60.0 + ra_s / 3600.0;
                    let ra_rad = ra_hours * 15.0 * (std::f64::consts::PI / 180.0);
                    let dec_d_str = parts[5];
                    let sign = if dec_d_str.starts_with('-') {
                        -1.0
                    } else {
                        1.0
                    };
                    let dec_d: f64 = dec_d_str.parse()?;
                    let dec_m: f64 = parts[6].parse()?;
                    let dec_s: f64 = parts[7].parse()?;
                    let dec_deg = sign * (dec_d.abs() + dec_m / 60.0 + dec_s / 3600.0);
                    let dec_rad = dec_deg.to_radians();
                    sources.push(Source {
                        name1,
                        name2,
                        ra_rad,
                        dec_rad,
                    });
                }
            }
            ParseSection::Sked => {
                let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
                if parts.contains(&"PREOB") && parts.contains(&"MIDOB") && parts.contains(&"POSTOB")
                {
                    if let Some(start_pos) = parts
                        .iter()
                        .position(|s| s.len() == 11 && s.chars().all(char::is_numeric))
                    {
                        let source_name = parts[0].to_string();
                        let start_str = parts[start_pos];
                        let duration_sec: i64 = parts[start_pos + 1].parse()?;
                        let year: i32 = 2000 + start_str[0..2].parse::<i32>()?;
                        let day_of_year: u32 = start_str[2..5].parse()?;
                        let hour: u32 = start_str[5..7].parse()?;
                        let minute: u32 = start_str[7..9].parse()?;
                        let second: u32 = start_str[9..11].parse()?;
                        if let Some(date) = NaiveDate::from_yo_opt(year, day_of_year) {
                            if let Some(datetime) = date.and_hms_opt(hour, minute, second) {
                                schedule.push(Observation {
                                    source_name,
                                    start_time: datetime,
                                    duration_sec,
                                });
                            }
                        }
                    }
                }
            }
            ParseSection::None => {}
        }
    }
    Ok(DrgData { sources, schedule })
}

fn get_default_stations() -> Vec<Station> {
    vec![
        Station {
            name: "KASHIM34".to_string(),
            pos: [-3997650.05799, 3276690.07124, 3724278.43114],
        },
        Station {
            name: "HITACH32".to_string(),
            pos: [-3961788.9740, 3243597.4920, 3790597.6920],
        },
        Station {
            name: "TAKAHA32".to_string(),
            pos: [-3961881.8250, 3243372.4800, 3790687.4490],
        },
        Station {
            name: "YAMAGU32".to_string(),
            pos: [-3502544.587, 3950966.235, 3566381.192],
        },
        Station {
            name: "YAMAGU34".to_string(),
            pos: [-3502567.576, 3950885.734, 3566449.115],
        },
    ]
}

// --- Main Execution ---
fn main() -> Result<(), eframe::Error> {
    let cli = Cli::parse();
    let all_stations = get_default_stations();

    let selected_station = match all_stations.iter().find(|s| s.name == cli.station) {
        Some(s) => s.clone(),
        None => {
            eprintln!(
                "Error: Station '{}' not found in the internal list.",
                cli.station
            );
            std::process::exit(1);
        }
    };

    let drg_data = match parse_drg_file(&cli.drg_file) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error parsing DRG file: {}", e);
            std::process::exit(1);
        }
    };

    if drg_data.schedule.is_empty() {
        eprintln!("Error: No schedule found in DRG file.");
        std::process::exit(1);
    }

    let t0 = drg_data
        .schedule
        .iter()
        .map(|obs| obs.start_time)
        .min()
        .unwrap();
    let max_time = drg_data
        .schedule
        .iter()
        .map(|obs| obs.start_time + Duration::seconds(obs.duration_sec))
        .max()
        .unwrap();

    if let Some(output_path) = cli.output.as_deref() {
        match write_full_track_file(
            output_path,
            &selected_station,
            &drg_data.sources,
            &drg_data.schedule,
            t0,
            max_time,
        ) {
            Ok(row_count) => eprintln!("Wrote {} full-track rows to {}", row_count, output_path),
            Err(e) => {
                eprintln!("Error writing {}: {}", output_path, e);
                std::process::exit(1);
            }
        }
    }

    if cli.terminal {
        print_terminal_observation_data(&selected_station, &drg_data.sources, &drg_data.schedule);
        return Ok(());
    }

    let min_x = 0.0;
    let max_x = (max_time - t0).num_seconds() as f64 / 3600.0;

    let range = max_x - min_x;
    let margin = range * 0.05; // 5% margin
    let x_axis_bounds = [min_x - margin, max_x + margin];

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1280.0, 720.0]),
        ..Default::default()
    };

    eframe::run_native(
        "DRG Uptime Plot",
        options,
        Box::new(move |cc| {
            let mut style = (*cc.egui_ctx.style()).clone();
            for (_text_style, font_id) in style.text_styles.iter_mut() {
                font_id.size *= 1.5;
            }
            cc.egui_ctx.set_style(style);
            Ok(Box::new(DrgPlotApp::new(
                selected_station,
                drg_data,
                x_axis_bounds,
                t0,
                max_time,
            )))
        }),
    )
}
