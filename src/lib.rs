//! World Generation tool which implements the `Generator trait` defined by `Robotic-Lib`
//!
//! # Features
//!
//! - Seed based generation.
//!     - All random elements of the world generation follow the provided seed. This means that every world is 100% reproducible.
//! - Elevation generation.
//!     - The world is generated starting from an elevation map which defines the shape of the terrain.
//! - Biomes generation.
//!     - Biomes are generated based on altitude and temperature.
//!     - Biomes list:
//!         - DeepWater
//!         - ShallowWater
//!         - Beach
//!         - Desert
//!           - In deserts you can fine lava lakes and fire patches.
//!         - Plain
//!         - Forest
//!         - Hill
//!         - Mountain
//!         - SnowyMountain
//! - Rivers generation
//!     - Rivers generation is based on altitude. Rivers spawn in the mountains and follow the altitude to end in a lake or the sea.
//! - Streets generation.
//!     - Streets spawn in the world connecting random points. This gives the robot a street infrastracture to use and start from.
//! - Telepors generation.
//!     - Teleports spawn in random locations.
//! - Content generation
//!     - Each content has a different spawn probability and for some contents the probability also changes based on the biome.
//!     - Contents:
//!         - `Tree` spawns in `Forest`, `Hill` and `Mountain`. Less trees spawn as the elevation increase
//!         - `Rock` spawns in `Plain`, `Hill`, `Mountain` and `SnowyMountain`. More rocks spawn as the elevation increase.
//!         - `Fish` spawns in `ShallowWater` and `DeepWater`. More fish spawn as the water depth increase.
//!     - `Garbage`, `Coins`, `Bins`, `Crates`, `Markets` spawn randomly in the world with different probabilities
//! - Weather forecast generation
//!     - The weather forecast is generated choosing a random weather for each day.
//!     - Weather types: `Sunny`, `Foggy`, `Rainy`, `TrentinoSnow`, `TropicalMonsoon`
//!     - `always_sunny` is a parameter that allow to generate an always sunny weather.
//! - Spawnpoint generation
//!     - The spawn point is generated randomly but with biome preference to make the start easier.
//!     - In order of preference: `Plain`, `Beach`, `Forest`
//!
//! # Example
//!
//! ```
//! use robotics_lib::world::world_generator::Generator;
//! use midgard::world_generator::WorldGenerator;
//! use midgard::world_generator::WorldGeneratorParameters;
//! use midgard::world_visualizer::WorldVisualizer;
//!
//! fn main() {
//!     // Define the WorldGenerator parameters using the dedicated struct
//!     let params = WorldGeneratorParameters {
//!         seed: 15, // fixed seed
//!         world_size: 200, // smaller world
//!         amount_of_rivers: None, // disable rivers
//!         amount_of_streets: Some(1.2), // more streets
//!         ..Default::default() // the rest of the parameters keep their default value
//!     };
//!
//!     // Instantiate the WorldGenerator with the parameters
//!     let mut world_generator = WorldGenerator::new(params);
//!     // Generate the world
//!     let (world, _spawn_point, _weather, _max_score, _score_table) = world_generator.gen();
//!
//!     // Use 'WorldVisualizer::visualize' to render the world at the specified resolution
//!     WorldVisualizer::visualize(world, 600);
//! }
//! ```

/// World Generator
pub mod world_generator;

/// A simple 2D visualizer to render the generated world
pub mod world_visualizer;
