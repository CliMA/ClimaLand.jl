# Sixel

![Julia version](https://img.shields.io/badge/julia-%3E%3D%201.6-blue)
[![Build Status](https://github.com/johnnychen94/Sixel.jl/workflows/CI/badge.svg)](https://github.com/johnnychen94/Sixel.jl/actions)
[![Coverage](https://codecov.io/gh/johnnychen94/Sixel.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/johnnychen94/Sixel.jl)

Encode the image into a [sixel][sixel_format_wiki] control sequence and vice versa. If your
terminal supports this format, then you can get a nice visualization of it.

This package, although itself complete, is probably not the most convinient way to use for normal
users. Package authors that aims to support sixel for various image-like types (e.g., image, video,
gif, graphic, plot, latex) are the potential targeted users of this package.

The core functionality of this package is supported by the C library [libsixel].

## Functions

This package exports two functions: `sixel_encode` and `sixel_decode`.

- `sixel_encode` converts the input array into sixel format sequence.
- `sixel_decode` converts the input sixel format sequence into colorant array.

## Terminals that support sixel

One important thing about sixel is that not all terminals support sixel control sequence. The
following is an incomplete list of terminals that support sixel.

- macOS: [iTerm2], [contour], [mlterm] and [wezterm]
- Linux: [mlterm], [contour], [xterm], [konsole] and [wezterm]
- Windows: [windows_terminal], [contour], [mintty], [msys2], [mlterm] and [wezterm]

> Above these I only manually test [iTerm2], [mlterm] and [wezterm].

> To run `xterm` with sixel support, use `$ xterm -xrm "XTerm*decTerminalID: vt340" -xrm "XTerm*numColorRegisters: 256"`.

A demo in iTerm2:

<img src="https://user-images.githubusercontent.com/8684355/118017901-b98f4d80-b389-11eb-93cf-0f99c169121d.png" alt="demo" width="1024" height="800">


It is the client terminal that renders the image. Hence if you're using any sixel-enabled terminal (e.g., iTerm2) and connect to other servers via ssh, you
still get sixel supported. Sixel is nothing special but a string sequence.

Unfortunately, there are some famous widely used advanced terminal/emulator that do not support
sixel (yet):

- [tmux does not support sixel](https://github.com/tmux/tmux/issues/1613#issuecomment-559940608).
- [kitty] does not support sixel; it has its own image protocol.

For more information, you can also read the [Terminal requirements](https://github.com/libsixel/libsixel#terminal-requirements) section in the [libsixel] repo.

Sixel.jl provides a function to test if your terminal supports it: `Sixel.is_sixel_supported()`. If
your terminal actually supports sixel and it returns `false`, please open an issue for it.

<!-- URLs -->

[windows_terminal]: https://github.com/microsoft/terminal
[contour]: https://github.com/christianparpart/contour
[iTerm2]: https://iterm2.com/
[kitty]: https://sw.kovidgoyal.net/kitty/
[libsixel]: https://github.com/libsixel/libsixel
[mlterm]: https://sourceforge.net/projects/mlterm/
[mintty]: https://github.com/mintty/mintty
[msys2]: https://www.msys2.org/
[sixel_format_wiki]: https://en.wikipedia.org/wiki/Sixel
[tmux]: https://github.com/tmux/tmux
[xterm]: https://invisible-island.net/xterm/
[konsole]: https://konsole.kde.org/
[wezterm]: https://wezfurlong.org/wezterm/
