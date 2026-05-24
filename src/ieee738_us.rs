/// Represents the physical, material, and electrical properties of a conductor.
#[derive(Debug, Clone, Copy)]
pub struct Conductor {
    /// Outer diameter of the conductor (ft)
    pub diameter: f64,
    /// Absorptivity of the conductor (0.0 to 1.0)
    pub absorptivity: f64,
    /// Emissivity of the conductor (0.0 to 1.0)
    pub emissivity: f64,
    /// Low Temperature (C) for resistance interpolation
    pub t_low: f64,
    /// High Temperature (C) for resistance interpolation
    pub t_high: f64,
    /// Resistance at Low Temperature (Ohms / ft)
    pub r_low: f64,
    /// Resistance at High Temperature (Ohms / ft)
    pub r_high: f64,
}

impl Conductor {
    /// Returns a standard Drake ACSR (795 kcmil 26/7) conductor profile.
    /// Uses standard diameter of 1.108 inches (0.092333 ft),
    /// resistance of 0.0220833 Ohms/kft (2.20833e-5 Ohms/ft) at 25°C,
    /// and 0.0263258 Ohms/kft (2.63258e-5 Ohms/ft) at 75°C.
    pub fn drake_acsr() -> Self {
        Self {
            diameter: 1.108 / 12.0,
            absorptivity: 0.8,
            emissivity: 0.8,
            t_low: 25.0,
            t_high: 75.0,
            r_low: 2.20833e-5,
            r_high: 2.63258e-5,
        }
    }

    /// Returns a standard Hawk ACSR (477 kcmil 26/7) conductor profile.
    /// Uses standard diameter of 0.858 inches (0.0715 ft),
    /// resistance of 0.036 Ohms/kft (3.6e-5 Ohms/ft) at 25°C,
    /// and 0.044 Ohms/kft (4.4e-5 Ohms/ft) at 75°C.
    pub fn hawk_acsr() -> Self {
        Self {
            diameter: 0.858 / 12.0,
            absorptivity: 0.8,
            emissivity: 0.8,
            t_low: 25.0,
            t_high: 75.0,
            r_low: 3.6e-5,
            r_high: 4.4e-5,
        }
    }
}

/// Represents the ambient weather and environmental conditions.
#[derive(Debug, Clone, Copy)]
pub struct Environment {
    /// Ambient temperature (C)
    pub ambient_temperature: f64,
    /// Wind speed (ft/s)
    pub wind_speed: f64,
    /// Wind angle (Degrees, 0 to 90)
    pub wind_angle_deg: f64,
    /// Height of conductor above sea level (ft)
    pub elevation: f64,
}

/// Represents the solar positioning and radiation inputs.
#[derive(Debug, Clone, Copy)]
pub struct Solar {
    /// Solar radiation (W/ft^2). If < 0.0, the solar radiation is calculated.
    pub solar_radiation: f64,
    /// Month of the year (1-12)
    pub month: i32,
    /// Day of the month (1-31)
    pub day_of_month: i32,
    /// Hour of the day (0.0-24.0, e.g. 11.0 is 11:00 AM)
    pub hour_of_day: f64,
    /// Latitude (Decimal Degrees)
    pub latitude_deg: f64,
    /// Line azimuth (Decimal Degrees, e.g. 90.0 for E-W line)
    pub line_azimuth_deg: f64,
    /// Clear atmosphere (true) vs industrial (false)
    pub atmosphere_clear: bool,
}


/// Returns convective_heat_loss Qc (Watts / ft)
/// # Arguments

/// * `ambient_temperature` - T_a: Degrees (C)
/// * `wind_speed` - V_w: Wind Speed (ft/s)
/// * `wind_angle_deg` - Wind Angle (Degrees) 0 to 90
/// * `elevation` - H_e: Height of conductor above sea level (ft)
/// * `conductor_temperature` - T_s: Conductor Surface Temperature (C)
/// * `diameter` - D_0: Outer diameter of the conductor (ft)
pub fn convective_heat_loss(
    ambient_temperature: f64,
    wind_speed: f64,
    wind_angle_deg: f64,
    elevation: f64,
    conductor_temperature: f64,
    diameter: f64,
) -> f64 {
    let pi = std::f64::consts::PI;

    // Limit to within 0-90.
    let wind_angle_deg_limited = 90.0 - (wind_angle_deg % 180.0 - 90.0).abs();
    let wind_angle_rad = wind_angle_deg_limited * (pi / 180.0);

    // Equation 6, Tfilm W/ft (degrees C)
    let tfilm = (conductor_temperature + ambient_temperature) / 2.0;

    // Absolute Viscosity of Air (m_f), (lb/ft*h)
    // dynamic_viscosity
    // Equation 13b
    let uf = 0.00353 * (tfilm + 273.15).powf(1.5) / (tfilm + 383.4);

    // air_density (lb/ft^3)
    // Equation 14b
    let pf = (0.080695 - 2.901e-6 * elevation + 3.7e-11 * elevation.powi(2)) / (1.0 + 0.00367 * tfilm);

    // Equation 4a Section 4.4.3.1, page 11.
    let kangle = 1.194
        - wind_angle_rad.cos()
        + 0.194 * (2.0 * wind_angle_rad).cos()
        + 0.368 * (2.0 * wind_angle_rad).sin();

    // Equation 2c
    let nre = diameter
        * pf
        * (wind_speed * 60.0 * 60.0) // Because dynamic_viscosity is in lb/ft-hr, we must convert wind speed to ft/hr.
        / uf;

    // thermal_conductivity_of_air
    // Equation 15b
    let kf = 7.388e-3 + 2.279e-5 * tfilm - 1.343e-9 * tfilm.powi(2);

    // Section 4.4.3.2, eq 5a 5b, page 12
    // qc0 = natural_convection
    let qc0 = 1.825
        * pf.powf(0.5)
        * diameter.powf(0.75)
        * (conductor_temperature - ambient_temperature).powf(1.25);

    // Equation 3a
    let qc1 = kangle
        * (1.01 + 1.35 * nre.powf(0.52))
        * kf
        * (conductor_temperature - ambient_temperature);

    // Equation 3b
    let qc2 = kangle
        * 0.754 * nre.powf(0.6)
        * kf
        * (conductor_temperature - ambient_temperature);

    // IEEE 738 recommends taking max of 3a / 3b results.
    // The convective heat loss is the bigger of forced and natural convection
    // From section 4.4.3 in the standard, page 10.
    f64::max(qc0, f64::max(qc1, qc2))
}

/// Returns radiated_heat_loss Qr (Watts / ft)
/// # Arguments
/// * `ambient_temperature` - T_a: Degrees (C)
/// * `conductor_temperature` - T_s: Conductor Surface Temperature (C)
/// * `emissivity` - ε: Epsilon, Emissivity of conductor (0.0 to 1.0)
/// * `diameter` - D_0: Outer diameter of the conductor (ft)
pub fn radiated_heat_loss(
    ambient_temperature: f64,
    conductor_temperature: f64,
    emissivity: f64,
    diameter: f64,
) -> f64 {
    // Section 4.4.4, eq 7a 7b, page 12
    1.656
        * diameter
        * emissivity
        * (
            ((conductor_temperature + 273.0) / 100.0).powi(4)
            - ((ambient_temperature + 273.0) / 100.0).powi(4)
        )
}

/// Calculates day of year. 
/// # Arguments
/// * `month` - Month January (1) to December (12)
/// * `day_of_month` - Day of Month, 1 to 31
pub fn day_of_year(month: i32, day_of_month: i32) -> i32 {
    let days_in_month = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    
    let mut result = day_of_month;
    for i in 1..month {
        result += days_in_month[i as usize];
    }
    result
}

/// Returns solar_heat_gain Qs (Watts / ft)
/// # Arguments
/// * `solar_radiation` - w/ft^2, or <0 if it should be calculated via month/day/hour
/// * `month` - 1 (January) to 12 (December)
/// * `day_of_month` - Day of Month (1-31)
/// * `hour_of_day` - Hour of Day, 0 to 23 (e.g. 11:00 AM => 11)
/// * `latitude_deg` - Lat: Latitude (Decimal Degrees)
/// * `line_azimuth_deg` - Z_l: If line runs E-W => 90 Degrees
/// * `elevation` - H_e: Height of conductor above sea level (ft)
/// * `atmosphere_clear` - Clear? (True) Industrial? (False)
/// * `absorptivity` - α: Alpha, Absorptivity of conductor (0.0 to 1.0)
/// * `diameter` - D_0: Outer diameter of the conductor (ft)
pub fn solar_heat_gain(
    solar_radiation: f64,
    month: i32,
    day_of_month: i32,
    hour_of_day: f64,
    latitude_deg: f64, 
    line_azimuth_deg: f64,
    elevation: f64, 
    atmosphere_clear: bool,
    absorptivity: f64,
    diameter: f64, 
) -> f64 {
    // If solar radiation is already specified, immediately return the value.
    if solar_radiation >= 0.0 {
        return absorptivity * solar_radiation * diameter;
    }

    // Constants
    let pi = std::f64::consts::PI;

    let day_of_year = day_of_year(month, day_of_month);

    let latitude_rad = latitude_deg * (pi / 180.0);

    // Hour angle relative to noon, 15*(Time-12), at 11AM, Time = 11 and the Hour angle= –15 deg 
    let w_deg = (hour_of_day - 12.0) * 15.0;
    let w_rad = w_deg * (pi / 180.0);

    // Table 3 - Atmosphere condition coefficients
    let (a, b, c, d, e, f, g) = match atmosphere_clear {
        true => (-3.9241, 5.9276, -1.7856e-1, 3.223e-3, -3.3549e-5, 1.8053e-7, -3.7868e-10),
        false => (4.9408, 1.3208, 6.1444e-2, -2.9411e-3, 5.07752e-5, -4.03627e-7, 1.22967e-9),
    };

    // Table H.5 - Solar heat multiplying factors, Ksolar for high altitudes
    let mult = match elevation {
        _ if elevation > 15000.0 => 1.3,
        _ if elevation > 10000.0 => 1.25,
        _ if elevation > 5000.0 => 1.15,
        _ => 1.0,
    };

    // Equation 16b - 23.4583 more precisely from Annex A
    let p_rad = (((284.0 + (day_of_year as f64)) / 365.0) * 360.0) * (pi / 180.0);
    let delta_rad = (23.4583 * p_rad.sin()) * (pi / 180.0);

    // Equation 16a
    let hc_rad = (latitude_rad.cos() * delta_rad.cos() * w_rad.cos() + latitude_rad.sin() * delta_rad.sin()).asin();
    // Limit to 0-90 range. Convert to degrees.
    let hc_deg = hc_rad * (180.0 / pi);

    // Equation 18 - Total solar and sky radiated heat intensity
    let qs = a + b * hc_deg + c * hc_deg.powi(2) + d * hc_deg.powi(3) + e * hc_deg.powi(4) + f * hc_deg.powi(5) + g * hc_deg.powi(6);

    // Equation 20 - Solar altitude correction factor
    let ksolar = 1.0 + 3.5e-5 * elevation - 1.0e-9 * elevation.powi(2);

    // Equation 8 - Total solar and sky radiated heat intensity corrected for elevation
    // Qs sometimes can compute as less than 0, if the sun is down. The lowest heating you can have is 0.
    let qse = f64::max(qs,0.0) * mult * ksolar;

    // Equation 17b
    let x = w_rad.sin() / ((latitude_rad.sin() * w_rad.cos() - latitude_rad.cos() * delta_rad.tan()));

    let cc_deg = 
        if -180.0 <= w_deg && w_deg < 0.0 {
            if x >= 0.0 { 0.0 } 
            else { 180.0 }
        } else {
            if x < 0.0 { 180.0 } 
            else { 360.0 }
        };

    let cc_rad = cc_deg * (pi / 180.0);

    // Azimuth of line
    let zl_rad = line_azimuth_deg * (pi / 180.0);

    // Azimuth of sun
    let zc_rad = cc_rad + (x).atan();

    // Equation 9 - Effective angle of incidence of the sun’s rays
    let theta = (hc_rad.cos() * (zc_rad - zl_rad).cos()).acos();

    // Compute solar_heat_flux
    absorptivity * qse * (theta).sin() * diameter
}

/// Returns resistance, adjusted to given conductor_temperature.
/// # Arguments
/// * `conductor_temperature` - T_s: Conductor Surface Temperature (C)
/// * `t_low` - Low Temperature, Degrees C
/// * `t_high` - High Temperature, Degrees C
/// * `r_low` - Resistance at Low Temperature, Ohms / ft
/// * `r_high` - Resistance at High Temperature, Ohms / ft
pub fn adjust_r(
    conductor_temperature: f64, 
    t_low: f64, 
    t_high: f64, 
    r_low: f64, 
    r_high: f64, 
) -> f64 {
    // Equation 10
    let ohms_per_c: f64 = (r_high - r_low) / (t_high - t_low);
    (ohms_per_c * (conductor_temperature - t_low)) + r_low
}

/// Returns thermal_rating (Amps)
/// # Arguments
/// * `conductor` - The Conductor properties parameter object
/// * `env` - The Environment conditions parameter object
/// * `solar` - The Solar positioning/radiation parameter object
/// * `conductor_temperature` - T_s: Conductor Surface Temperature (C)
pub fn thermal_rating(
    conductor: &Conductor,
    env: &Environment,
    solar: &Solar,
    conductor_temperature: f64,
) -> f64 {

    if conductor_temperature < env.ambient_temperature {
        return 0.0;
    }

    let qc = convective_heat_loss(
        env.ambient_temperature,
        env.wind_speed,
        env.wind_angle_deg,
        env.elevation,
        conductor_temperature,
        conductor.diameter,
    );

    let qr = radiated_heat_loss(
        env.ambient_temperature,
        conductor_temperature,
        conductor.emissivity,
        conductor.diameter,
    );

    let qs: f64 = solar_heat_gain(
        solar.solar_radiation,
        solar.month,
        solar.day_of_month,
        solar.hour_of_day,
        solar.latitude_deg,
        solar.line_azimuth_deg,
        env.elevation,
        solar.atmosphere_clear,
        conductor.absorptivity,
        conductor.diameter,
    );

    let r = adjust_r(
        conductor_temperature,
        conductor.t_low,
        conductor.t_high,
        conductor.r_low,
        conductor.r_high,
    );

    if qc + qr - qs < 0.0 {
        // The ambient temperature + solar heating, has brought the conductor to a higher temperature than the specified MOT "conductor_temperature"
        return 0.0;
    }

    ((qc + qr - qs) / r).powf(0.5)
}

/// Returns calculated_temperature (C) based on input conditions
/// # Arguments
/// * `conductor` - The Conductor properties parameter object
/// * `env` - The Environment conditions parameter object
/// * `solar` - The Solar positioning/radiation parameter object
/// * `current` - Current (amps)
/// * `tolerance` - Tolerance on result (amps)
pub fn calculated_temperature(
    conductor: &Conductor,
    env: &Environment,
    solar: &Solar,
    current: f64,
    tolerance: f64,
) -> f64 {
    if current < 0.0 {
        return 0.0;
    }

    let mut lower_bound: f64 = env.ambient_temperature;
    let mut upper_bound: f64 = 256.0;
    let target_y: f64 = current;

    // Increase upper_bound until y(upper_bound) exceeds target_y or it becomes very large
    loop {
        let thermal_rating_retval = thermal_rating(
            conductor,
            env,
            solar,
            upper_bound,
        );

        if thermal_rating_retval < target_y && upper_bound < f64::MAX / 2.0 {
            upper_bound *= 2.0;
        } else {
            break;
        }
    }

    // Bisection search with known upper_bound and lower_bound
    while upper_bound - lower_bound > tolerance {
        let mid = (lower_bound + upper_bound) / 2.0;
        let mid_y = thermal_rating(
            conductor,
            env,
            solar,
            mid,
        );

        if mid_y <= target_y {
            lower_bound = mid;
        } else {
            upper_bound = mid;
        }
    }

    // Return the midpoint of the final range
    (lower_bound + upper_bound) / 2.0
}

/// Returns conductor_temperature_rise (C)
/// # Arguments
/// * `conductor` - The Conductor properties parameter object
/// * `env` - The Environment conditions parameter object
/// * `solar` - The Solar positioning/radiation parameter object
/// * `conductor_temperature` - Initial Conductor Surface Temperature (C)
/// * `current` - Current (amps)
/// * `time_step` - Timestep (seconds)
/// * `steps` - Number of time steps to apply
/// * `heat_capacity` - m*Cp: Total heat capacity of conductor (J/(ft-°C))
pub fn conductor_temperature_rise(
    conductor: &Conductor,
    env: &Environment,
    solar: &Solar,
    conductor_temperature: f64,
    current: f64,
    time_step: f64,
    steps: i32,
    heat_capacity: f64,
) -> f64 {

    if conductor_temperature < env.ambient_temperature {
        return 0.0;
    }

    let mut final_temperature = conductor_temperature;

    for _ in 0..steps {
        let qc = convective_heat_loss(
            env.ambient_temperature,
            env.wind_speed,
            env.wind_angle_deg,
            env.elevation,
            final_temperature,
            conductor.diameter,
        );
        let qr = radiated_heat_loss(
            env.ambient_temperature,
            final_temperature,
            conductor.emissivity,
            conductor.diameter,
        );
        let qs: f64 = solar_heat_gain(
            solar.solar_radiation,
            solar.month,
            solar.day_of_month,
            solar.hour_of_day,
            solar.latitude_deg,
            solar.line_azimuth_deg,
            env.elevation,
            solar.atmosphere_clear,
            conductor.absorptivity,
            conductor.diameter,
        );
        let r = adjust_r(
            final_temperature,
            conductor.t_low,
            conductor.t_high,
            conductor.r_low,
            conductor.r_high,
        );
        let delta_t = ((r * current.powf(2.0)) + qs - qc - qr) * time_step / heat_capacity;
        final_temperature += delta_t;
    }

    final_temperature - conductor_temperature
}

/// Returns transient_rating (Amps)
/// # Arguments
/// * `conductor` - The Conductor properties parameter object
/// * `env` - The Environment conditions parameter object
/// * `solar` - The Solar positioning/radiation parameter object
/// * `conductor_temperature` - Initial Conductor Surface Temperature (C)
/// * `conductor_temperature_max` - Max Final Conductor Surface Temperature (C)
/// * `time_step` - Timestep (seconds)
/// * `steps` - Number of time steps to apply
/// * `tolerance` - Tolerance on result (amps)
/// * `heat_capacity` - m*Cp: Total heat capacity of conductor (J/(ft-°C))
pub fn transient_rating(
    conductor: &Conductor,
    env: &Environment,
    solar: &Solar,
    conductor_temperature: f64,
    conductor_temperature_max: f64,
    time_step: f64,
    steps: i32,
    tolerance: f64,
    heat_capacity: f64,
) -> f64 {

    if conductor_temperature_max < conductor_temperature {
        return 0.0;
    }

    // Assume the rating is somewhere between 0 to 4096A.
    let mut lower_bound: f64 = 0.0;
    let mut upper_bound: f64 = 4096.0;
    // delta_t_max
    let target_y: f64 = conductor_temperature_max - conductor_temperature; 

    // Increase upper_bound until y(upper_bound) exceeds target_y or it becomes very large
    while conductor_temperature_rise(
        conductor,
        env,
        solar,
        conductor_temperature,
        upper_bound,
        time_step,
        steps,
        heat_capacity,
    ) < target_y && upper_bound < f64::MAX / 2.0 {
        upper_bound *= 2.0;
    }

    // Bisection search with known upper_bound and lower_bound
    while upper_bound - lower_bound > tolerance {
        let mid = (lower_bound + upper_bound) / 2.0;
        let mid_y = conductor_temperature_rise(
            conductor,
            env,
            solar,
            conductor_temperature,
            mid,
            time_step,
            steps,
            heat_capacity,
        );

        if mid_y < target_y {
            lower_bound = mid;
        } else {
            upper_bound = mid;
        }
    }

    // Return the midpoint of the final range
    (lower_bound + upper_bound) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_drake_acsr_thermal_rating() {
        let conductor = Conductor::drake_acsr();
        let env = Environment {
            ambient_temperature: 40.0,
            wind_speed: 2.0,
            wind_angle_deg: 90.0,
            elevation: 0.0,
        };
        let solar = Solar {
            solar_radiation: -1.0, // calculated
            month: 6,
            day_of_month: 10,
            hour_of_day: 11.0,
            latitude_deg: 30.0,
            line_azimuth_deg: 90.0,
            atmosphere_clear: true,
        };

        // Calculate thermal rating at max operating temperature (MOT) of 100°C
        let rating = thermal_rating(&conductor, &env, &solar, 100.0);
        
        // Under these conditions (standard solar, 40°C ambient, 2 ft/s crosswind, MOT 100°C),
        // the Drake rating should be positive and roughly around 800 - 1100 Amps.
        assert!(rating > 0.0);
        assert!(rating > 800.0 && rating < 1100.0, "Drake rating: {}", rating);
        
        // Verify calculated_temperature returns 100°C given the thermal rating
        let temp = calculated_temperature(&conductor, &env, &solar, rating, 0.01);
        assert!((temp - 100.0).abs() < 0.1);
    }

    #[test]
    fn test_hawk_acsr_thermal_rating() {
        let conductor = Conductor::hawk_acsr();
        let env = Environment {
            ambient_temperature: 25.0,
            wind_speed: 2.0,
            wind_angle_deg: 90.0,
            elevation: 0.0,
        };
        let solar = Solar {
            solar_radiation: 1000.0, // explicit solar radiation
            month: 1,
            day_of_month: 1,
            hour_of_day: 12.0,
            latitude_deg: 0.0,
            line_azimuth_deg: 0.0,
            atmosphere_clear: true,
        };

        // Hawk has smaller diameter and higher resistance than Drake.
        // Therefore, it should have a lower thermal rating under similar conditions.
        let drake = Conductor::drake_acsr();
        
        let hawk_rating = thermal_rating(&conductor, &env, &solar, 75.0);
        let drake_rating = thermal_rating(&drake, &env, &solar, 75.0);
        
        assert!(hawk_rating > 0.0);
        assert!(drake_rating > hawk_rating, "Drake rating ({}) should be higher than Hawk rating ({}) due to larger diameter & lower resistance", drake_rating, hawk_rating);
    }
}
