fn main() {
    match dockq::cli::run() {
        Ok(response) => {
            if let Some((path, contents)) = response.json_file {
                if let Err(error) = std::fs::write(&path, contents) {
                    eprintln!(
                        "Error: could not write JSON results to {}: {error}",
                        path.display()
                    );
                    std::process::exit(1);
                }
            }
            print!("{}", response.stdout);
        }
        Err(error) => {
            eprintln!("Error: {error}");
            std::process::exit(1);
        }
    }
}
