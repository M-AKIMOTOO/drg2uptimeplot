use astro::{coords, ecliptic, sun, time};
use chrono::{DateTime, Datelike, Duration, NaiveDate, NaiveDateTime, TimeZone, Timelike, Utc};
use clap::Parser;
use eframe::egui;
use egui::{Color32, Stroke};
use egui_plot::{Corner, GridMark, Line, LineStyle, Plot, PlotBounds, PlotPoints, Points};
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

    /// Initial station for AZ/EL calculations (also selectable in the GUI)
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

struct StationPolarTracks {
    name: String,
    sky_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
    drg_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
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

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum AppTab {
    UptimePlot,
    UptimePlot2,
    PolarPlot1,
    PolarPlot2,
    YI,
}

// --- Plotting App ---
struct DrgPlotApp {
    selected_station: Station,
    available_stations: Vec<Station>,
    sources: Vec<Source>,
    schedule: Vec<Observation>,
    sky_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
    drg_ut_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
    plot_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)>,
    color_map: HashMap<String, Color32>,
    x_axis_bounds: [f64; 2],
    t0: NaiveDateTime,
    t_end: NaiveDateTime,
    selected_tab: AppTab,
    reset_plot_bounds: bool,
    hidden_targets: HashMap<AppTab, HashSet<String>>,
    yi_yamaguchi32: StationPolarTracks,
    yi_yamaguchi34: StationPolarTracks,
    yi_yamaguchi34_offset_enu: [f64; 2],
}

impl DrgPlotApp {
    fn new(
        station: Station,
        available_stations: Vec<Station>,
        yi_stations: [Station; 2],
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

        let yi_yamaguchi32 = calculate_station_polar_tracks(
            &yi_stations[0],
            &drg_data.sources,
            &drg_data.schedule,
            t0,
            t_end,
        );
        let yi_yamaguchi34 = calculate_station_polar_tracks(
            &yi_stations[1],
            &drg_data.sources,
            &drg_data.schedule,
            t0,
            t_end,
        );
        let yi_yamaguchi34_offset_enu = station_offset_enu(&yi_stations[0], &yi_stations[1]);

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
            selected_station: station,
            available_stations,
            sources: drg_data.sources,
            schedule: drg_data.schedule,
            sky_segments,
            drg_ut_segments,
            plot_segments,
            color_map,
            x_axis_bounds,
            t0,
            t_end,
            selected_tab: AppTab::UptimePlot,
            reset_plot_bounds: false,
            hidden_targets: HashMap::new(),
            yi_yamaguchi32,
            yi_yamaguchi34,
            yi_yamaguchi34_offset_enu,
        }
    }
}

impl eframe::App for DrgPlotApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let current_station_name = self.selected_station.name.clone();
        let mut selected_station_name = current_station_name.clone();
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.selected_tab, AppTab::UptimePlot, "UptimePlot1");
                ui.selectable_value(&mut self.selected_tab, AppTab::UptimePlot2, "UptimePlot2");
                ui.selectable_value(&mut self.selected_tab, AppTab::PolarPlot1, "PolarPlot1");
                ui.selectable_value(&mut self.selected_tab, AppTab::PolarPlot2, "PolarPlot2");
                ui.selectable_value(&mut self.selected_tab, AppTab::YI, "YI");

                ui.separator();
                ui.label("Station:");
                egui::ComboBox::from_id_salt("station_selector")
                    .selected_text(&selected_station_name)
                    .show_ui(ui, |ui| {
                        for station in &self.available_stations {
                            ui.selectable_value(
                                &mut selected_station_name,
                                station.name.clone(),
                                &station.name,
                            );
                        }
                    });

                if ui.button("Reset Zoom").clicked() {
                    self.reset_plot_bounds = true;
                }
            });
        });

        if selected_station_name != current_station_name {
            if let Some(station) = self
                .available_stations
                .iter()
                .find(|station| station.name == selected_station_name)
                .cloned()
            {
                self.selected_station = station;
                self.recalculate_selected_station_data();
                self.reset_plot_bounds = true;
            }
        }

        egui::SidePanel::right("target_legend_panel")
            .default_width(190.0)
            .resizable(true)
            .show(ctx, |ui| self.ui_target_legend(ui));

        egui::CentralPanel::default().show(ctx, |ui| match self.selected_tab {
            AppTab::UptimePlot => self.ui_uptime_plot_tab(ui),
            AppTab::UptimePlot2 => self.ui_uptime_plot2_tab(ui),
            AppTab::PolarPlot1 => self.ui_polar_plot1_tab(ui),
            AppTab::PolarPlot2 => self.ui_polar_plot2_tab(ui),
            AppTab::YI => self.ui_yi_tab(ui),
        });
        self.reset_plot_bounds = false;
    }
}

impl DrgPlotApp {
    fn recalculate_selected_station_data(&mut self) {
        let sky_segments = calculate_ut_sky_segments(
            &self.selected_station,
            &self.sources,
            &self.schedule,
            self.t0,
            self.t_end,
        );
        let drg_ut_segments = calculate_observation_segments_ut(
            &self.selected_station,
            &self.sources,
            &self.schedule,
        );
        let mut plot_segments = calculate_observation_segments(
            &self.selected_station,
            &self.sources,
            &self.schedule,
            self.t0,
        );
        plot_segments.push(calculate_sun_segments(
            &self.selected_station,
            self.t0,
            self.t_end,
            false,
        ));

        self.sky_segments = sky_segments;
        self.drg_ut_segments = drg_ut_segments;
        self.plot_segments = plot_segments;
    }

    fn target_legend_entries(&self) -> Vec<(String, Color32)> {
        let mut names = Vec::new();
        match &self.selected_tab {
            AppTab::UptimePlot | AppTab::PolarPlot1 => {
                names.extend(self.plot_segments.iter().map(|(name, _, _)| name.clone()));
            }
            AppTab::UptimePlot2 | AppTab::PolarPlot2 => {
                names.extend(self.sky_segments.iter().map(|(name, _, _)| name.clone()));
                names.extend(
                    self.drg_ut_segments
                        .iter()
                        .filter(|(name, _, _)| name != "Sun")
                        .map(|(name, _, _)| name.clone()),
                );
            }
            AppTab::YI => {
                for tracks in [&self.yi_yamaguchi32, &self.yi_yamaguchi34] {
                    names.extend(tracks.sky_segments.iter().map(|(name, _, _)| name.clone()));
                    names.extend(tracks.drg_segments.iter().map(|(name, _, _)| name.clone()));
                }
            }
        }
        names.sort();
        names.dedup();
        names
            .into_iter()
            .filter_map(|name| self.color_map.get(&name).map(|color| (name, *color)))
            .collect()
    }

    fn ui_target_legend(&mut self, ui: &mut egui::Ui) {
        ui.heading("Targets");
        if matches!(
            self.selected_tab,
            AppTab::UptimePlot2 | AppTab::PolarPlot2 | AppTab::YI
        ) {
            ui.small("Dotted: sky   Solid: DRG");
        }
        ui.separator();

        let entries = self.target_legend_entries();
        let hidden_targets = self.hidden_targets.entry(self.selected_tab).or_default();
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                for (name, color) in entries {
                    let mut visible = !hidden_targets.contains(&name);
                    ui.horizontal(|ui| {
                        if ui.checkbox(&mut visible, "").changed() {
                            if visible {
                                hidden_targets.remove(&name);
                            } else {
                                hidden_targets.insert(name.clone());
                            }
                        }
                        let (swatch, _) =
                            ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                        ui.painter().rect_filled(swatch, 2.0, color);
                        ui.label(name);
                    });
                }
            });
    }

    fn target_is_hidden(&self, name: &str) -> bool {
        self.hidden_targets
            .get(&self.selected_tab)
            .is_some_and(|targets| targets.contains(name))
    }

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
            );

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
            );

        plot_az.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max(
                    [self.x_axis_bounds[0], 0.0],
                    [self.x_axis_bounds[1], 360.0],
                ));
            }
            for (name, az_points, _) in &self.plot_segments {
                if self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    for points in split_finite_segments(az_points) {
                        plot_ui
                            .line(Line::new(name.clone(), PlotPoints::from(points)).color(*color));
                    }
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
                if self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    for points in split_finite_segments(el_points) {
                        plot_ui
                            .line(Line::new(name.clone(), PlotPoints::from(points)).color(*color));
                    }
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
            );

        let plot_el = Plot::new("el_plot_ut24")
            .width(ui.available_width())
            .default_x_bounds(0.0, 24.0)
            .height(ui.available_height() / 2.0)
            .y_axis_label("Elevation (deg)")
            .y_axis_min_width(70.0)
            .allow_drag(true)
            .allow_zoom(true)
            .allow_scroll(true)
            .include_y(0.0)
            .include_y(90.0)
            .set_margin_fraction(egui::Vec2::ZERO)
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
            );

        plot_az.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([0.0, 0.0], [24.0, 360.0]));
            }

            for (name, az_points, _) in &self.sky_segments {
                if self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    for points in split_finite_segments(az_points) {
                        plot_ui.line(
                            Line::new(format!("{} (sky)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(1.0)
                                .style(LineStyle::dotted_dense()),
                        );
                    }
                }
            }

            for (name, az_points, _) in &self.drg_ut_segments {
                if name == "Sun" || self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    for points in split_finite_segments(az_points) {
                        plot_ui.line(
                            Line::new(format!("{} (DRG)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(3.0),
                        );
                    }
                }
            }
        });

        ui.add_space(-10.0);

        plot_el.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([0.0, 0.0], [24.0, 90.0]));
            }

            for (name, _, el_points) in &self.sky_segments {
                if self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let visible_points = filter_elevation_points(el_points, 5.0);
                    for points in split_finite_segments(&visible_points) {
                        plot_ui.line(
                            Line::new(format!("{} (sky)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(1.0)
                                .style(LineStyle::dotted_dense()),
                        );
                    }
                }
            }

            for (name, _, el_points) in &self.drg_ut_segments {
                if name == "Sun" || self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let visible_points = filter_elevation_points(el_points, 5.0);
                    for points in split_finite_segments(&visible_points) {
                        plot_ui.line(
                            Line::new(format!("{} (DRG)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(3.0),
                        );
                    }
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
            .allow_scroll(true); // Added for interactivity

        plot.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([-1.0, -1.0], [1.0, 1.0]));
            }
            draw_polar_grid(plot_ui, 0.0);

            for (name, az_points, el_points) in &self.plot_segments {
                if self.target_is_hidden(name) {
                    continue;
                }
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

            plot_ui.points(
                Points::new(
                    "Antenna (AZ=244, EL=20)",
                    PlotPoints::from(vec![azel_to_polar_point(244.0, 20.0)]),
                )
                .color(Color32::WHITE)
                .radius(6.0),
            );
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
            .allow_scroll(true);

        plot.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([-1.0, -1.0], [1.0, 1.0]));
            }
            draw_polar_grid(plot_ui, 0.0);

            for (name, az_points, el_points) in &self.sky_segments {
                if self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let points = azel_to_polar_points(az_points, el_points);
                    for points in split_finite_segments(&points) {
                        plot_ui.line(
                            Line::new(format!("{} (sky)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(1.0)
                                .style(LineStyle::dotted_dense()),
                        );
                    }
                }
            }

            for (name, az_points, el_points) in &self.drg_ut_segments {
                if name == "Sun" || self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let points = azel_to_polar_points(az_points, el_points);
                    for points in split_finite_segments(&points) {
                        plot_ui.line(
                            Line::new(format!("{} (DRG)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(3.0),
                        );
                    }
                }
            }

            plot_ui.points(
                Points::new(
                    "Antenna (AZ=244, EL=20)",
                    PlotPoints::from(vec![azel_to_polar_point(244.0, 20.0)]),
                )
                .color(Color32::WHITE)
                .radius(6.0),
            );
        });
    }

    fn ui_yi_tab(&self, ui: &mut egui::Ui) {
        let area = ui.available_rect_before_wrap();
        ui.allocate_rect(area, egui::Sense::hover());

        let [east, north] = self.yi_yamaguchi34_offset_enu;
        let distance = east.hypot(north);
        let bearing = east.atan2(north).to_degrees();
        ui.painter().text(
            area.left_top() + egui::vec2(8.0, 8.0),
            egui::Align2::LEFT_TOP,
            format!(
                "YAMAGU34 relative to YAMAGU32: E={east:.1} m, N={north:.1} m, {distance:.1} m at {bearing:.1} deg"
            ),
            egui::FontId::proportional(14.0),
            Color32::GRAY,
        );

        let chart_area = area.shrink2(egui::vec2(8.0, 28.0));
        let chart_size = (chart_area.width() * 0.42)
            .min(chart_area.height() * 0.42)
            .max(1.0);
        let baseline_direction = egui::vec2(east as f32, -(north as f32)).normalized();
        let chart_offset = baseline_direction * (chart_size * 1.52);
        let chart_size = egui::vec2(chart_size, chart_size);
        let yamagu32_rect =
            egui::Rect::from_center_size(chart_area.center() - chart_offset * 0.5, chart_size);
        let yamagu34_rect =
            egui::Rect::from_center_size(chart_area.center() + chart_offset * 0.5, chart_size);

        self.draw_yi_station_polar(ui, yamagu32_rect, "yi_yamaguchi32", &self.yi_yamaguchi32);
        self.draw_yi_station_polar(ui, yamagu34_rect, "yi_yamaguchi34", &self.yi_yamaguchi34);
    }

    fn draw_yi_station_polar(
        &self,
        parent_ui: &mut egui::Ui,
        rect: egui::Rect,
        plot_id: &'static str,
        tracks: &StationPolarTracks,
    ) {
        let mut station_ui = parent_ui.new_child(egui::UiBuilder::new().max_rect(rect));
        station_ui.set_clip_rect(rect);
        station_ui.label(&tracks.name);

        let plot = Plot::new(plot_id)
            .width(station_ui.available_width())
            .height(station_ui.available_height())
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
            .allow_scroll(true);

        plot.show(&mut station_ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([-1.0, -1.0], [1.0, 1.0]));
            }
            draw_polar_grid(plot_ui, 0.0);

            for (name, az_points, el_points) in &tracks.sky_segments {
                if self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let points = azel_to_polar_points(az_points, el_points);
                    for points in split_finite_segments(&points) {
                        plot_ui.line(
                            Line::new(format!("{} (sky)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(1.0)
                                .style(LineStyle::dotted_dense()),
                        );
                    }
                }
            }

            for (name, az_points, el_points) in &tracks.drg_segments {
                if name == "Sun" || self.target_is_hidden(name) {
                    continue;
                }
                if let Some(color) = self.color_map.get(name) {
                    let points = azel_to_polar_points(az_points, el_points);
                    for points in split_finite_segments(&points) {
                        plot_ui.line(
                            Line::new(format!("{} (DRG)", name), PlotPoints::from(points))
                                .color(*color)
                                .width(3.0),
                        );
                    }
                }
            }

            plot_ui.points(
                Points::new(
                    "Antenna (AZ=244, EL=20)",
                    PlotPoints::from(vec![azel_to_polar_point(244.0, 20.0)]),
                )
                .color(Color32::WHITE)
                .radius(6.0),
            );
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
            if !azimuth.is_finite() || !elevation.is_finite() || elevation < 5.0 {
                return [f64::NAN, f64::NAN];
            }
            azel_to_polar_point(azimuth, elevation)
        })
        .collect()
}

fn azel_to_polar_point(azimuth: f64, elevation: f64) -> [f64; 2] {
    let angle = (90.0 - azimuth).to_radians();
    let radius = (90.0 - elevation) / 90.0;
    [radius * angle.cos(), radius * angle.sin()]
}

fn filter_elevation_points(points: &[[f64; 2]], minimum_elevation: f64) -> Vec<[f64; 2]> {
    points
        .iter()
        .map(|point| {
            if point[1].is_finite() && point[1] >= minimum_elevation {
                *point
            } else {
                [point[0], f64::NAN]
            }
        })
        .collect()
}

fn split_finite_segments(points: &[[f64; 2]]) -> Vec<Vec<[f64; 2]>> {
    let mut segments = Vec::new();
    let mut current = Vec::new();

    for &point in points {
        if point[0].is_finite() && point[1].is_finite() {
            current.push(point);
        } else if !current.is_empty() {
            segments.push(std::mem::take(&mut current));
        }
    }

    if !current.is_empty() {
        segments.push(current);
    }
    segments
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

fn calculate_station_polar_tracks(
    station: &Station,
    sources: &[Source],
    schedule: &[Observation],
    t0: NaiveDateTime,
    t_end: NaiveDateTime,
) -> StationPolarTracks {
    StationPolarTracks {
        name: station.name.clone(),
        sky_segments: calculate_sky_segments(station, sources, schedule, t0, t_end),
        drg_segments: calculate_observation_segments(station, sources, schedule, t0),
    }
}

fn station_offset_enu(origin: &Station, target: &Station) -> [f64; 2] {
    let origin_ecef = ECEF::new(origin.pos[0], origin.pos[1], origin.pos[2]);
    let origin_wgs84: WGS84<f64> = origin_ecef.into();
    let longitude = origin_wgs84.longitude_radians();
    let latitude = origin_wgs84.latitude_radians();
    let [dx, dy, dz] = [
        target.pos[0] - origin.pos[0],
        target.pos[1] - origin.pos[1],
        target.pos[2] - origin.pos[2],
    ];

    let east = -longitude.sin() * dx + longitude.cos() * dy;
    let north = -latitude.sin() * longitude.cos() * dx - latitude.sin() * longitude.sin() * dy
        + latitude.cos() * dz;
    [east, north]
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
    let yi_stations = [
        all_stations
            .iter()
            .find(|station| station.name == "YAMAGU32")
            .expect("YAMAGU32 must be present in the internal station list")
            .clone(),
        all_stations
            .iter()
            .find(|station| station.name == "YAMAGU34")
            .expect("YAMAGU34 must be present in the internal station list")
            .clone(),
    ];
    let available_stations = all_stations.clone();

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
                available_stations,
                yi_stations,
                drg_data,
                x_axis_bounds,
                t0,
                max_time,
            )))
        }),
    )
}

#[cfg(test)]
mod tests {
    use super::{get_default_stations, station_offset_enu};

    #[test]
    fn yamagu34_is_northeast_of_yamagu32() {
        let stations = get_default_stations();
        let yamagu32 = stations
            .iter()
            .find(|station| station.name == "YAMAGU32")
            .unwrap();
        let yamagu34 = stations
            .iter()
            .find(|station| station.name == "YAMAGU34")
            .unwrap();
        let [east, north] = station_offset_enu(yamagu32, yamagu34);

        assert!(east > 0.0);
        assert!(north > 0.0);
        assert!((east - 70.6).abs() < 0.2);
        assert!((north - 81.5).abs() < 0.2);
    }
}
