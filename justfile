pkg  := "shazam"
desc := "I/O for GMRT ring buffers, with the power of SHAZAM!"

alias d := docs
alias c := clean
alias i := install
alias u := uninstall

# List available commands.
default:
    #!/usr/bin/env python
    from rich.table import Table
    from rich.panel import Panel
    from rich.console import Console

    console = Console()

    grid = Table.grid(expand=True, padding=(0, 2, 0, 2))
    grid.add_column(justify="left", style="bold")
    grid.add_column(justify="right", style="italic")

    grid.add_row("clean ([i]c[/i])", "Clean up")
    grid.add_row("install ([i]i[/i])", "Install")
    grid.add_row("docs ([i]d[/i])", "Build docs.")
    grid.add_row("uninstall ([i]u[/i])", "Uninstall")

    console.print(
        Panel(
            grid,
            padding=2,
            expand=False,
            title="[b]{{pkg}}[/b]: [i]{{desc}}[/i]",
        )
    )

# Clean up.
@clean:
    echo "Cleaning..."
    rm -rf docs/build/*
    fd -I tmp -x rm -rf
    fd -I dist -x rm -rf
    fd -I .eggs -x rm -rf
    fd -I .cache -x rm -rf
    fd -I -e pyc -x rm -rf
    fd -I .coverage -x rm -rf
    fd -I .mypy_cache -x rm -rf
    fd -I __pycache__ -x rm -rf
    fd -I .pytest_cache -x rm -rf

# Install.
@install: && clean
    echo "Installing..."
    pip install --no-build-isolation -Ceditable.rebuild=true -ve .

# Uninstall.
@uninstall: && clean
    echo "Uninstalling {{pkg}}..."
    pip uninstall {{pkg}}
    rm -rf python/{{pkg}}.egg-info
    rm -rf python/{{pkg}}/_version.py

# Build docs.
@docs:
    echo "Building docs for {{pkg}}..."
    sphinx-reload docs
