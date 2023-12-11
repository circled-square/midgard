//! World Generator
//! 
//! World Generation tool wich implements the `Generator trait` defined by `Robotic-Lib`
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
//!         - `Trees` spawn in `Forest`, `Hill` and `Mountain`. Less trees spawn as the elevation increase
//!         - `Rocks` spawn in `Plain`, `Hill`, `Mountain` and `SnowyMountain`. More rocks spawn as the elevation increase.
//!         - `Fish` spawn in `ShallowWater` and `DeepWater`. More fish spawn as the water depth increase.
//!     - `Garbage`, `Coins`, `Bins`, `Crates`, `Markets` spawn randomly in the world with different probabilities
//! - Weather forecast generation
//!     - The weather forecast is generated choosing a random weather for each day.
//!     - Weather types: `Sunny`, `Foggy`, `Rainy`, `TrentinoSnow`, `TropicalMonsoon`
//!     - `always_sunny` is a parameter that allow to generate an always sunny weather.
//! - Spawnpoint generation
//!     - The spawn point is generated randomly but with biome preference to make the start easier. 
//!     - In order of preference: `Plain`, `Beach`, `Forest`
//! 
//! # Examples
//! 
//! ```
//! // Import the 'Generator' trait and  'WorldGeneratorParameters' to use the 'WorldGenerator'
//! use midgard::world_generator::WorldGeneratorParameters;
//! use robotics_lib::world::world_generator::Generator;
//! 
//! // Import the "WorldGenerator"
//! use midgard::world_generator::WorldGenerator;
//! 
//! // Import the visualizer if you want to view a 2D render of your world
//! use midgard::world_visualizer::WorldVisualizer;
//! 
//! # fn main() {
//! // Define the World Generator parameters using the dedicated struct
//! let params = WorldGeneratorParameters {
//!     world_size: 300,
//!     ..Default::default()
//! };
//! 
//! // Instantiate the World Generator with the static method 'new', passing the parameters
//! let mut world_generator = WorldGenerator::new(params);
//! let (world, (_spawn_x, _spawn_y), _weather, _max_score, _score_table) = world_generator.gen();
//! 
//! // Use the 'visualize' method to render the generated world
//! // the 2nd parameter is the window resolution and the 3rd is the scaling
//! WorldVisualizer::visualize(world, 600, 2);
//! # }
//! ```

/// World Generator
pub mod world_generator;

/// A simple 2D visualizer to render the generated world
pub mod world_visualizer;
