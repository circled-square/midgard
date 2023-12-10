use robotics_lib::world::tile::Content;
use std::collections::HashMap;

pub struct WorldGeneratorParameters {
    pub seed: u64,
    pub world_size: usize,
    pub always_sunny: bool,
    pub time_progression_minutes: u8,
    pub starting_hour: u8,
    pub world_scale: f64,
    pub amount_of_rivers: Option<f64>,
    pub amount_of_streets: Option<f64>,
    pub amount_of_teleports: Option<f64>,
    pub score_table: Option<HashMap<Content, f32>>,
    pub max_score: f32,
    pub contents_radii: ContentsRadii,
}
impl Default for WorldGeneratorParameters {
    fn default() -> Self {
        Self {
            seed: rand::random(),
            world_size: 300,
            always_sunny: false,
            time_progression_minutes: 10,
            starting_hour: 8,
            world_scale: 180.0,
            amount_of_rivers: Some(1.0),
            amount_of_streets: Some(1.0),
            amount_of_teleports: Some(1.0),
            score_table: None,
            max_score: 1000.0,
            contents_radii: ContentsRadii::default(),
        }
    }
}

pub struct ContentsRadii {
    pub trees_in_forest: u64,
    pub trees_in_hill: u64,
    pub trees_in_mountain: u64,
    pub rocks_in_plains: u64,
    pub rocks_in_hill: u64,
    pub rocks_in_mountain: u64,
    pub fish_in_shallow_water: u64,
    pub fish_in_deep_water: u64,
    pub garbage: u64,
    pub coins: u64,
    pub garbage_bins: u64,
    pub crates: u64,
    pub markets: u64,
}
impl Default for ContentsRadii {
    fn default() -> Self {
        Self {
            trees_in_forest: 3,
            trees_in_hill: 4,
            trees_in_mountain: 5,
            rocks_in_plains: 5,
            rocks_in_hill: 4,
            rocks_in_mountain: 3,
            fish_in_shallow_water: 5,
            fish_in_deep_water: 4,
            garbage: 10,
            coins: 20,
            garbage_bins: 20,
            crates: 40,
            markets: 50,
        }
    }
}
