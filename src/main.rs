use robotics_lib::world::worldgenerator::Generator;
use robotics_lib::world::tile::Tile;
use robotics_lib::world::environmental_conditions::EnvironmentalConditions;
use robotics_lib::world::environmental_conditions::WeatherType;
use robotics_lib::world::tile::TileType;
use robotics_lib::world::tile::Content;

struct WorldGenerator {}
impl Generator for WorldGenerator {
    fn gen(&mut self) -> (Vec<Vec<Tile>>, (usize, usize), EnvironmentalConditions, f32) {
        let world = vec![vec![Tile {
            tile_type: TileType::Grass,
            content: Content::Tree(1),
            elevation: 0,
        }; 10]; 10];
        let pos = (0,0);
        let env_conditions = EnvironmentalConditions::new(&[WeatherType::Sunny], 30, 8);
        let max_score = 100.;

        (world, pos, env_conditions, max_score)
    }
}

fn main() {
    println!("Hello, world!");
}
