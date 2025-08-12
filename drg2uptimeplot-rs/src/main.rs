mod utils;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use chrono::{Duration, NaiveDate, NaiveDateTime, TimeZone, Utc};
use clap::Parser;
use eframe::egui;
use egui::{Color32, Stroke};
use egui_plot::{GridMark, Legend, Line, Plot, PlotPoints, Corner, PlotBounds, Points};

// --- Command Line Arguments Definition ---
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Path to the DRG schedule file
    drg_file: String,

    /// Name of the station to use for calculations
    #[arg(long, default_value = "YAMAGU32")]
    station: String,
}

// --- Data Structures ---
#[derive(Debug, Clone)]
struct Source {
    name: String,
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
    PolarPlot,
}

// --- Plotting App ---
struct DrgPlotApp {
    plot_segments: Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>) >,
    color_map: HashMap<String, Color32>,
    x_axis_bounds: [f64; 2],
    t0: NaiveDateTime,
    selected_tab: AppTab,
    reset_plot_bounds: bool,
}

impl DrgPlotApp {
    fn new(station: Station, drg_data: DrgData, x_axis_bounds: [f64; 2], t0: NaiveDateTime) -> Self {
        let plot_segments = calculate_observation_segments(&station, &drg_data.sources, &drg_data.schedule, t0);
        
        let mut color_map = HashMap::new();
        let palette = [
            Color32::from_rgb(100, 143, 255), // Blue
            Color32::from_rgb(120, 255, 120), // Green
            Color32::from_rgb(255, 100, 100), // Red
            Color32::from_rgb(255, 180, 80),  // Orange
            Color32::from_rgb(240, 120, 240), // Magenta
            Color32::from_rgb(130, 255, 255), // Cyan
        ];

        let mut unique_sources = drg_data.sources.iter().map(|s| s.name.clone()).collect::<Vec<_>>();
        unique_sources.sort();
        unique_sources.dedup();

        for (i, name) in unique_sources.iter().enumerate() {
            color_map.insert(name.clone(), palette[i % palette.len()]);
        }

        Self {
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
                ui.selectable_value(&mut self.selected_tab, AppTab::UptimePlot, "Uptime Plot");
                ui.selectable_value(&mut self.selected_tab, AppTab::PolarPlot, "Polar Plot");

                ui.separator();

                if ui.button("Reset Zoom").clicked() {
                    self.reset_plot_bounds = true;
                }
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            match self.selected_tab {
                AppTab::UptimePlot => self.ui_uptime_plot_tab(ui),
                AppTab::PolarPlot => self.ui_polar_plot_tab(ui),
            }
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

        let az_pointer_formatter = |x: f64, y: f64| format!("Time: {:02}:{:02}\nAz: {:.1}°", x as u32, (x.fract() * 60.0) as u32, y);
        let el_pointer_formatter = |x: f64, y: f64| format!("Time: {:02}:{:02}\nEl: {:.1}°", x as u32, (x.fract() * 60.0) as u32, y);

        let plot_az = Plot::new("az_plot").width(ui.available_width()).height(ui.available_height() / 2.0)
            .y_axis_label("Azimuth (deg)")
            .y_axis_min_width(70.0) // Changed from 0.0 to 70.0 for uniform width
            .allow_drag(true).allow_zoom(true).allow_scroll(true)
            .include_y(0.0).include_y(360.0)
            .y_grid_spacer(|_input| (0..=12).map(|v| GridMark { value: (v * 30) as f64, step_size: 30.0 }).collect()) // Re-added
            .show_x(true) // Re-added, Temporarily set to true for testing
            .x_axis_label("") // Re-added
            .x_axis_formatter(|_, _| "".to_string()) // Re-added
            .y_axis_formatter(|m, _| format!("{:>3}", m.value as i32)).show_y(true) // Re-added, changed to 3-digit padded
            .coordinates_formatter(Corner::LeftTop, egui_plot::CoordinatesFormatter::new(move |p, _| az_pointer_formatter(p.x, p.y)))
            .legend(Legend::default());

        let plot_el = Plot::new("el_plot").width(ui.available_width()).height(ui.available_height() / 2.0)
            .y_axis_label("Elevation (deg)")
            .y_axis_min_width(70.0) // Changed from 67.0 to 70.0 for uniform width
            .allow_drag(true).allow_zoom(true).allow_scroll(true)
            .include_y(0.0).include_y(90.0)
            .x_axis_formatter(move |m, _| pointer_time_formatter(m.value)) // Re-added
            .y_axis_formatter(|m, _| format!("{:>2}", m.value as i32)).show_y(true) // Re-added, changed to 3-digit padded
            .y_grid_spacer(|_input| (0..=9).map(|v| GridMark { value: (v * 10) as f64, step_size: 10.0 }).collect()) // Re-added
            .coordinates_formatter(Corner::LeftTop, egui_plot::CoordinatesFormatter::new(move |p, _| el_pointer_formatter(p.x, p.y)))
            .legend(Legend::default());

        plot_az.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([self.x_axis_bounds[0], 0.0], [self.x_axis_bounds[1], 360.0]));
            }
            for (name, az_points, _) in &self.plot_segments {
                if let Some(color) = self.color_map.get(name) {
                    plot_ui.line(Line::new(name.clone(), PlotPoints::from(az_points.clone())).color(*color));
                }
            }
        });

        ui.add_space(-10.0);

        plot_el.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([self.x_axis_bounds[0], 0.0], [self.x_axis_bounds[1], 90.0]));
            }
            for (name, _, el_points) in &self.plot_segments {
                 if let Some(color) = self.color_map.get(name) {
                    plot_ui.line(Line::new(name.clone(), PlotPoints::from(el_points.clone())).color(*color));
                }
            }
        });
    }

    fn ui_polar_plot_tab(&mut self, ui: &mut egui::Ui) {
        let plot = Plot::new("polar_plot")
            .width(ui.available_width()).height(ui.available_height())
            .data_aspect(1.0).view_aspect(1.0) 
            .include_x(-1.0).include_x(1.0)
            .include_y(-1.0).include_y(1.0)
            .center_x_axis(true).center_y_axis(true)
            .show_x(false).show_y(false)
            .x_grid_spacer(|_input| vec![]).y_grid_spacer(|_input| vec![])
            .allow_drag(true).allow_zoom(true).allow_scroll(true) // Added for interactivity
            .legend(Legend::default());

        plot.show(ui, |plot_ui| {
            if self.reset_plot_bounds {
                plot_ui.set_plot_bounds(PlotBounds::from_min_max([-1.0, -1.0], [1.0, 1.0]));
            }
            for el_level in [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0] {
                let radius = (90.0 - el_level) / 90.0;
                if radius >= 0.0 {
                    let circle_points: PlotPoints = (0..=100).map(|i| {
                        let angle = i as f64 * 2.0 * std::f64::consts::PI / 100.0;
                        [radius * angle.cos(), radius * angle.sin()]
                    }).collect();
                    plot_ui.line(Line::new("", circle_points).stroke(Stroke::new(1.0, Color32::DARK_GRAY)));
                    if el_level != 90.0 {
                        let label_text = format!("{:.0}°", el_level);
                        let label_x = radius * (72.0f64).to_radians().cos();
                        let label_y = radius * (72.0f64).to_radians().sin();
                        plot_ui.text(egui_plot::Text::new("", egui_plot::PlotPoint::new(label_x, label_y), label_text).color(Color32::DARK_GRAY));
                    }
                }
            }

            for az_level in [0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0] {
                let angle_rad = (90.0f64 - az_level).to_radians();
                let line_points = vec![[0.0, 0.0], [angle_rad.cos(), angle_rad.sin()]] ;
                plot_ui.line(Line::new("", PlotPoints::from(line_points)).stroke(Stroke::new(1.0, Color32::DARK_GRAY)));
                let label_text = format!("{:.0}°", az_level);
                plot_ui.text(egui_plot::Text::new("", egui_plot::PlotPoint::new(angle_rad.cos() * 1.1, angle_rad.sin() * 1.1), label_text).color(Color32::DARK_GRAY));
            }

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
                        plot_ui.points(Points::new(name.clone(), PlotPoints::from(polar_points)).color(*color));
                    }
                }
            }
        });
    }
}

// --- Calculation Logic ---
fn calculate_observation_segments(station: &Station, sources: &[Source], schedule: &[Observation], t0: NaiveDateTime) -> Vec<(String, Vec<[f64; 2]>, Vec<[f64; 2]>)> {
    let mut new_plot_data = Vec::new();
    let ant_pos = station.pos;

    for obs in schedule {
        if let Some(source) = sources.iter().find(|s| s.name == obs.source_name) {
            let mut az_segment = Vec::new();
            let mut el_segment = Vec::new();

            let start_time = obs.start_time;
            let end_time = start_time + Duration::seconds(obs.duration_sec);

            let mut current_time = start_time;
            while current_time <= end_time {
                let duration_since_t0 = current_time.signed_duration_since(t0);
                let hour_float = duration_since_t0.num_seconds() as f64 / 3600.0;

                let datetime_utc = Utc.from_utc_datetime(&current_time);
                let (az, el, _) = utils::radec2azalt(ant_pos, datetime_utc, source.ra_rad, source.dec_rad);

                if el >= 0.0 {
                    az_segment.push([hour_float, az]);
                    el_segment.push([hour_float, el]);
                } else {
                    az_segment.push([hour_float, az]);
                    el_segment.push([hour_float, f64::NAN]);
                }
                current_time += Duration::minutes(1);
            }
            new_plot_data.push((source.name.clone(), az_segment, el_segment));
        }
    }
    new_plot_data
}

// --- File Parsing Logic ---
fn parse_drg_file<P: AsRef<Path>>(path: P) -> Result<DrgData, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sources = Vec::new();
    let mut schedule = Vec::new();
    enum ParseSection { None, Sources, Sked }
    let mut current_section = ParseSection::None;

    for line in reader.lines() {
        let line = line?;
        let trimmed_line = line.trim();
        if trimmed_line.starts_with('$') {
            current_section = match trimmed_line {
                "$SOURCES" => ParseSection::Sources,
                "$SKED" => ParseSection::Sked,
                _ => ParseSection::None
            };
            continue;
        }
        if trimmed_line.is_empty() || trimmed_line.starts_with('*') { continue; }

        match current_section {
            ParseSection::Sources => {
                let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
                if parts.len() >= 9 && parts[8] == "2000.0" {
                    let name = parts[0].to_string();
                    let ra_h: f64 = parts[2].parse()?;
                    let ra_m: f64 = parts[3].parse()?;
                    let ra_s: f64 = parts[4].parse()?;
                    let ra_hours = ra_h + ra_m / 60.0 + ra_s / 3600.0;
                    let ra_rad = ra_hours * 15.0 * (std::f64::consts::PI / 180.0);
                    let dec_d_str = parts[5];
                    let sign = if dec_d_str.starts_with('-') { -1.0 } else { 1.0 };
                    let dec_d: f64 = dec_d_str.parse()?;
                    let dec_m: f64 = parts[6].parse()?;
                    let dec_s: f64 = parts[7].parse()?;
                    let dec_deg = sign * (dec_d.abs() + dec_m / 60.0 + dec_s / 3600.0);
                    let dec_rad = dec_deg.to_radians();
                    sources.push(Source { name, ra_rad, dec_rad });
                }
            }
            ParseSection::Sked => {
                let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
                if parts.contains(&"PREOB") && parts.contains(&"MIDOB") && parts.contains(&"POSTOB") {
                    if let Some(start_pos) = parts.iter().position(|s| s.len() == 11 && s.chars().all(char::is_numeric)) {
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
                                 schedule.push(Observation { source_name, start_time: datetime, duration_sec });
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
        Station { name: "KASHIM34".to_string(), pos: [-3997650.05799, 3276690.07124, 3724278.43114] },
        Station { name: "HITACH32".to_string(), pos: [-3961788.9740, 3243597.4920, 3790597.6920] },
        Station { name: "TAKAHA32".to_string(), pos: [-3961881.8250, 3243372.4800, 3790687.4490] },
        Station { name: "YAMAGU32".to_string(), pos: [-3502544.587, 3950966.235, 3566381.192] },
        Station { name: "YAMAGU34".to_string(), pos: [-3502567.576, 3950885.734, 3566449.115] },
    ]
}

// --- Main Execution ---
fn main() -> Result<(), eframe::Error> {
    let cli = Cli::parse();
    let all_stations = get_default_stations();

    let selected_station = match all_stations.iter().find(|s| s.name == cli.station) {
        Some(s) => s.clone(),
        None => {
            eprintln!("Error: Station '{}' not found in the internal list.", cli.station);
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

    let t0 = drg_data.schedule.iter().map(|obs| obs.start_time).min().unwrap();
    let max_time = drg_data.schedule.iter().map(|obs| obs.start_time + Duration::seconds(obs.duration_sec)).max().unwrap();
    
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
            Ok(Box::new(DrgPlotApp::new(selected_station, drg_data, x_axis_bounds, t0)))
        }),
    )
}
