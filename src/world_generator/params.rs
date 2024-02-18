use robotics_lib::world::tile::Content;
use std::collections::HashMap;

/// Contains parameters passed to `world_generator::WorldGenerator` to tweak its behaviour
///
/// These parameters are used to personalize the world generation process, for example
/// setting the world size or its scaling, disabling certain features, or
/// changing the rarity of certain contents. For most use cases `WorldGeneratorParameters::default()`
/// should be ok, and it is recommended when setting parameters to start from a default instance.
///
/// # Examples
/// Users can simply use the default parameters:
/// ```
/// # use midgard::{*, params::*};
/// # use robotics_lib::world::world_generator::Generator;
/// let mut world_generator = WorldGenerator::new(WorldGeneratorParameters::default());
/// let (world, spawn_point, weather, max_score, score_table) = world_generator.gen();
/// ```
///
/// Or they can change them to their liking:
/// ```
/// # use midgard::{*, params::*};
/// # use robotics_lib::world::world_generator::Generator;
/// let params = WorldGeneratorParameters {
///     seed: 15, // fixed seed
///     world_size: 200, // smaller world
///     amount_of_rivers: None, // disable rivers
///     amount_of_streets: Some(1.2), // more streets
///     ..Default::default() // the rest of the parameters keep their default value
/// };
/// let mut world_generator = WorldGenerator::new(params);
/// let (world, spawn_point, weather, max_score, score_table) = world_generator.gen();
/// ```

#[derive(Clone)]
pub struct WorldGeneratorParameters {
    /// Seed used for world generation.
    pub seed: u64,

    /// length of the side of the world, which is always a square
    pub world_size: usize,

    /// if true disables weather generation and the weather will always be sunny
    pub always_sunny: bool,

    /// the number of days of weather that should be generated
    pub weather_forecast_length: u64,

    /// the amount of minutes that pass for each tick
    pub time_progression_minutes: u8,

    /// the starting hour, for example if `starting hour == 8` after generation the time will be 8:00 AM
    pub starting_hour: u8,

    /// the scaling of the world. a smaller scale will result in smaller mountains, valleys and lakes,
    /// and shorter distances between them
    pub world_scale: f64,


    /// Controls the amount of rivers generated.
    /// If set to `None` the river generation step will be skipped.
    pub amount_of_rivers: Option<f64>,

    /// Controls the amount of streets generated.
    /// If set to `None` the street generation step will be skipped.
    pub amount_of_streets: Option<f64>,

    /// Controls the amount of teleports generated.
    /// If set to `None` the teleport generation step will be skipped.
    pub amount_of_teleports: Option<f64>,

    /// Controls the maximum elevation, and scales all elevation accordingly.
    pub elevation_multiplier: Option<f64>,

    /// Controls the amount of each tile content to be spawned. See [ ContentsRadii ]
    pub contents_radii: ContentsRadii,

    /// Sets a custom score table.
    /// If set to `None` the default one provided by `robotics_lib` will be used.
    pub score_table: Option<HashMap<Content, f32>>,

    /// sets the maximum score the robot can earn
    pub max_score: f32,
}

impl Default for WorldGeneratorParameters {
    /// The default values are the following:
    /// ```
    /// # use midgard::{*, params::*};
    /// # WorldGeneratorParameters {
    /// seed: rand::random(),
    /// world_size: 300,
    /// always_sunny: false,
    /// weather_forecast_length: 7,
    /// time_progression_minutes: 10,
    /// starting_hour: 8,
    /// world_scale: 1.0,
    /// amount_of_rivers: Some(1.0),
    /// amount_of_streets: Some(1.0),
    /// amount_of_teleports: Some(1.0),
    /// elevation_multiplier: Some(4.0),
    /// contents_radii: ContentsRadii::default(),
    /// score_table: None,
    /// max_score: 1000.0,
    /// # };
    /// ```
    fn default() -> Self {
        Self {
            seed: rand::random(),
            world_size: 300,
            always_sunny: false,
            weather_forecast_length: 7,
            time_progression_minutes: 10,
            starting_hour: 8,
            world_scale: 1.0,
            amount_of_rivers: Some(1.0),
            amount_of_streets: Some(1.0),
            amount_of_teleports: Some(1.0),
            elevation_multiplier: Some(4.0),
            contents_radii: ContentsRadii::default(),
            score_table: None,
            max_score: 1000.0,
        }
    }
}


/// Controls the amount of each tile content to be spawned.
///
/// The values in the struct can be thought of as the "rarity" of each tile content (sometimes
/// specific to a particular biome); for example if `trees_in_forest == 3` and `trees_in_hill == 4`
/// that means that trees are more rare in hills than they are in forests.
///
/// What the numbers actually represent are the radii (or radiuses) of the Poisson distributions
/// used to generate the contents.
/// For most use cases `ContentsRadii::default()` should be ok, and it is recommended when setting
/// parameters to start from a default instance.
///
/// ***Note**: the user may notice that changing some of these values can cause a performance hit;
/// this is because the Poisson distributions used are cached to avoid generating the same more
/// than once. Setting them to a value different from all others in the struct trades this
/// performance benefit for some added customization. For small world sizes the performance hit
/// should be minimal, and it is up to the user to decide what to prioritize.*

#[derive(Clone)]
pub struct ContentsRadii {
    pub trees_in_forest: u64,
    pub trees_in_hill: u64,
    pub trees_in_mountain: u64,
    pub rocks_in_plains: u64,
    pub rocks_in_hill: u64,
    pub rocks_in_mountain: u64,
    pub bushes_in_plains: u64,
    pub fish_in_shallow_water: u64,
    pub fish_in_deep_water: u64,
    pub garbage: u64,
    pub coins: u64,
    pub garbage_bins: u64,
    pub crates: u64,
    pub markets: u64,
    pub banks: u64,
    pub buildings: u64,
    pub scarecrows: u64,
    pub jolly_blocks: u64
}

impl Default for ContentsRadii {
    /// The default values are the following:
    /// ```
    /// # midgard::params::ContentsRadii {
    /// trees_in_forest: 3,
    /// trees_in_hill: 4,
    /// trees_in_mountain: 5,
    /// rocks_in_plains: 5,
    /// rocks_in_hill: 4,
    /// rocks_in_mountain: 3,
    /// bushes_in_plains: 4,
    /// fish_in_shallow_water: 5,
    /// fish_in_deep_water: 4,
    /// garbage: 10,
    /// coins: 30,
    /// garbage_bins: 20,
    /// crates: 40,
    /// markets: 50,
    /// banks: 50,
    /// buildings: 50,
    /// scarecrows: 30,
    /// jolly_blocks: 50,
    /// # };
    /// ```
    fn default() -> Self {
        Self {
            trees_in_forest: 3,
            trees_in_hill: 4,
            trees_in_mountain: 5,
            rocks_in_plains: 5,
            rocks_in_hill: 4,
            rocks_in_mountain: 3,
            bushes_in_plains: 4,
            fish_in_shallow_water: 5,
            fish_in_deep_water: 4,
            garbage: 10,
            coins: 30,
            garbage_bins: 20,
            crates: 40,
            markets: 50,
            banks: 50,
            buildings: 50,
            scarecrows: 30,
            jolly_blocks: 50,
        }
    }
}
