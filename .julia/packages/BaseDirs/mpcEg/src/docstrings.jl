# Reference: https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html

Base.Docs.doc!(@__MODULE__, # Please let me know if there's an easier way.
               Base.Docs.Binding(@__MODULE__, Symbol(@__MODULE__)),
               Base.Docs.docstr(
                   """
**BaseDirs**

This module provides utilities to identify the appropriate locations for files.

The [`User`](@ref BaseDirs.User) and [`System`](BaseDirs.System) submodules
provide a number of accessor functions, which see.

There are also four combined accessor functions defined, namely: [`data`](@ref
BaseDirs.data), [`config`](@ref BaseDirs.config), [`fonts`](@ref
BaseDirs.fonts), and [`applications`](@ref BaseDirs.applications). These provide
a list of all relevant user and system directories.

As a convenience, the three user-specific accessors [`cache`](@ref
BaseDirs.User.cache), [`runtime`](@ref BaseDirs.User.runtime), and
[`state`](@ref BaseDirs.User.state), are available under `BaseDirs` as they have
no `System` equivalent they could be confused with.

Also see [`BaseDirs.Project`](@ref BaseDirs.Project) for information on how to
generate appropriate project paths.

!!! note
    This is essentially an implementation of the XDG (Cross-Desktop Group) directory
    specifications, with analogues for Windows and MacOS for cross-platform. More
    specifically, this is a hybrid of:
    - The [*XDG base directory*](https://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html)
      and the [*XDG user directory*](https://www.freedesktop.org/wiki/Software/xdg-user-dirs/) specifications on Linux
    - The [*Known Folder*](https://msdn.microsoft.com/en-us/library/windows/desktop/dd378457.aspx) API on Windows
    - The [*Standard Directories*](https://developer.apple.com/library/content/documentation/FileManagement/Conceptual/FileSystemProgrammingGuide/FileSystemOverview/FileSystemOverview.html#//apple_ref/doc/uid/TP40010672-CH2-SW6)
      guidelines on macOS

# Guidelines for appropriate directories

+ **Runtime** data is _volatile_, and only relevant to the current user session.
  It should be user-specific and non-essential. For example: sockets, named
  pipes, or lockfiles.
+ **Cache** data is non-essential and can be recreated _without data loss_. For
  example: thumbnails, compiled bytecode, font paths, or other cached data.
+ **State** data relates to application state, that should persist across _sessions_
  but doesn't need to be backed up. For example: logs, recently opened files, or
  other session data.
+ **Config** data is user-specific customisations to application behaviour. For
  example: theme settings, custom keybindings, media preferences, or other
  configuration.
+ **Data** content should be used for general application data. For example: saved
  game progress, library metadata, templates, or other
  user-specific/user-generated data.
+ **Bin** data is specifically executable files intended to be run by the user.

Essential information (that you'd want to include in backups) should be split
across **config** and **data**. Information in **state** may be nice to have,
but is non-essential.
""",
                   Dict(:path => joinpath(@__DIR__, @__FILE__),
                        :linenumber => @__LINE__,
                        :module => @__MODULE__)),
               Union{})

@doc """
     Project

A representation of a "Project", namely the essential components of naming
information used to produce platform-appropriate project paths.

    Project(name::Union{AbstractString, Module};
            org::AbstractString="julia", qualifier::AbstractString="lang")

The information needed, and the platforms that make use of it, are as follows:
- `name`, the name of the project (Linux, MacOS, Windows)
- `org` (`"julia"`), the organisation the project belongs to (MacOS, Windows)
- `qualifier` (`"org"`), the nature of the organisation, usually a TLD (MacOS)

The resulting "project path components" take one of the following forms:

| Platform | Project path form         |
|----------|---------------------------|
| Linux    | `"\$org/\$name"`            |
| MacOS    | `"\$qualifier.\$org.\$name"` |
| Windows  | `"\$org\\\$name"`           |
""" Project

# ---------

@doc """
**DATA_HOME** (`XDG_DATA_HOME`)

The single base directory relative to which user-specific data files should be
written.

**Default values**

| Linux            | MacOS                          | Windows          |
|------------------|--------------------------------|------------------|
| `~/.local/share` | `~/Library/ApplicationSupport` | `RoamingAppData` |
""" DATA_HOME

@doc """
**DATA_DIRS** (`XDG_DATA_DIRS`)

The set of *preference ordered* base directories relative to which data files
should be searched.

**Default values**

| Linux              | MacOS                         | Windows       |
|--------------------|-------------------------------|---------------|
| `/usr/local/share` | `/Library/ApplicationSupport` | `ProgramData` |
| `/usr/share`       |                               |               |
""" DATA_DIRS

@doc """
**CONFIG_HOME** (`XDG_CONFIG_HOME`)

The single base directory relative to which user-specific configuration files
should be written.

**Default values**

| Linux             | MacOS                          | Windows          |
|-------------------|--------------------------------|------------------|
| `~/.local/config` | `~/Library/ApplicationSupport` | `RoamingAppData` |
""" CONFIG_HOME

@doc """
**CONFIG_DIRS** (`XDG_CONFIG_DIRS`)

The set of *preference ordered* base directories relative to which data files
should be searched.

**Default values**

| Linux      | MacOS                         | Windows       |
|------------|-------------------------------|---------------|
| `/etc/xdg` | `/Library/ApplicationSupport` | `ProgramData` |
""" CONFIG_DIRS

@doc """
**STATE_HOME** (`XDG_STATE_HOME`)

The single base directory relative to which user-specific state data
should be written.

This should contain  state data that should persist between (application)
restarts, but that is not important or portable enough to the user that it
should be stored in `DATA_HOME`. It may contain:
- actions history (logs, history, recently used files, …)
- current state of the application that can be reused on a restart (view, layout, open files, undo history, …)

**Default values**

| Linux            | MacOS                          | Windows        |
|------------------|--------------------------------|----------------|
| `~/.local/state` | `~/Library/ApplicationSupport` | `LocalAppData` |
""" STATE_HOME

@doc """

**BIN_HOME** (`XDG_BIN_HOME`)

The single base directory relative to which user-specific executables should be
written.

**Default values**

| Linux          | MacOS\\*                         | Windows\\*             |
|----------------|--------------------------------|----------------------|
| `~/.local/bin` | `~/.local/bin`                 | `\\bin`               |
|                | `/opt/local/bin`               | `RoamingAppData\\bin` |
|                | `/usr/local/bin`               | `AppData\\bin`        |
|                |                                | current working dir  |

\\* The first of these directories that exists is used.

!!! warning
    This is not yet standardised by the XDG, see
    [freedesktop.org/xdg-specs#14](https://gitlab.freedesktop.org/xdg/xdg-specs/-/issues/14)
    for more information.
""" BIN_HOME

@doc """
**CACHE_HOME** (`XDG_CACHE_HOME`)

The single base directory relative to which user-specific non-essential (cached)
data should be written.

**Default values**

| Linux      | MacOS              | Windows              |
|------------|--------------------|----------------------|
| `~/.cache` | `~/Library/Caches` | `LocalAppData\\cache` |

""" CACHE_HOME

@doc """
**RUNTIME_DIR** (`XDG_RUNTIME_DIR`)

The single base directory relative to which user-specific runtime files and
other file objects should be placed. . Applications should use this directory
for communication and synchronization purposes and should not place larger files
in it.

**Default values**

| Linux            | MacOS                          | Windows        |
|------------------|--------------------------------|----------------|
| `/run/user/\$UID` | `~/Library/ApplicationSupport` | `LocalAppData` |
""" RUNTIME_DIR

# ---------

@doc """
**DESKTOP_DIR** (`XDG_DESKTOP_DIR`)

The user's desktop directory.
""" DESKTOP_DIR

@doc """
**DOWNLOAD_DIR** (`XDG_DOWNLOAD_DIR`)

The user's downloads directory.
""" DOWNLOAD_DIR

@doc """
**DOCUMENTS_DIR** (`XDG_DOCUMENTS_DIR`)

The user's documents directory.
""" DOCUMENTS_DIR

@doc """
**MUSIC_DIR** (`XDG_MUSIC_DIR`)

The user's music directory.
""" MUSIC_DIR

@doc """
**PICTURES_DIR** (`XDG_PICTURES_DIR`)

The user's pictures directory.
""" PICTURES_DIR

@doc """
**VIDEOS_DIR** (`XDG_VIDEOS_DIR`)

The user's videos directory.
""" VIDEOS_DIR

@doc """
**TEMPLATES_DIR** (`XDG_TEMPLATES_DIR`)

The user's templates directory.
""" TEMPLATES_DIR

@doc """
**PUBLICSHARE_DIR** (`XDG_PUBLICSHARE_DIR`)

The user's public directory.
""" PUBLICSHARE_DIR

@doc """
**APPLICATIONS_DIRS**

A list of locations in which application files may be found/placed.
""" APPLICATIONS_DIRS

@doc """
**FONTS_DIRS**

A list of locations in which font files may be found.
""" FONTS_DIRS

# ---------

@doc """
**BaseDirs.User**

This module containes accessor functions for user-specific directories.

### Base directory accessors

The functions [`data`](@ref BaseDirs.User.data), [`config`](@ref
BaseDirs.User.config), [`state`](@ref BaseDirs.User.state), [`cache`](@ref
BaseDirs.User.cache), and [`runtime`](@ref BaseDirs.User.runtime) produce a
`String`.

### User directory accessors

The functions [`desktop`](@ref BaseDirs.User.desktop), [`downloads`](@ref
BaseDirs.User.downloads), [`music`](@ref BaseDirs.User.music), [`videos`](@ref
BaseDirs.User.videos), [`templates`](@ref BaseDirs.User.templates),
[`public`](@ref BaseDirs.User.public) produce a `String`.

### Other accessors

The functions [`fonts`](@ref BaseDirs.User.fonts) and [`applications`](@ref
BaseDirs.User.applications) produce a `Vector{String}`.

!!! note
    Unlike the [`System`](@ref BaseDirs.System) and "combined" (`BaseDirs.*`)
    accessors, the *Base* and *User* accessors here return a single directory
    (`String`).
""" User

@doc """
$(Internals.accessordoc(:data, :DATA_HOME, name="user configuration"))
""" User.data

@doc """
$(Internals.accessordoc(:config, :CONFIG_HOME, name="user configuration"))
""" User.config

@doc """
$(Internals.accessordoc(:state, :STATE_HOME, name="state data"))
""" User.state

@doc """
$(Internals.accessordoc(:bin, :BIN_HOME, name="executables"))

!!! note "Special behaviour"
    When `create` is `true` and the path referrs to a file, `chmod` is called to
    ensure that all users who can read the file can execute it.
""" User.bin

@doc """
$(Internals.accessordoc(:cache, :CACHE_HOME, name="cached data"))
""" User.cache

@doc """
$(Internals.accessordoc(:state, :STATE_HOME, name="runtime information"))
""" User.runtime

@doc """
$(Internals.accessordoc(:fonts, name="user fonts", plural=true))
""" User.fonts

@doc """
$(Internals.accessordoc(:applications, name="user applications", plural=true))
""" User.applications

# ---------

@doc """
    desktop(parts...) -> String

Join the desktop directory with zero or more path components (`parts`).

The desktop directory is based on the variable `BaseDirs.DESKTOP_DIR`, which see.
""" User.desktop


@doc """
    downloads(parts...) -> String

Join the downloads directory with zero or more path components (`parts`).

The downloads directory is based on the variable `BaseDirs.DOWNLOADS_DIR`, which see.
""" User.downloads

@doc """
    documents(parts...) -> String

Join the documents directory with zero or more path components (`parts`).

The documents directory is based on the variable `BaseDirs.DOCUMENTS_DIR`, which see.
""" User.documents

@doc """
    music(parts...) -> String

Join the music directory with zero or more path components (`parts`).

The music directory is based on the variable `BaseDirs.MUSIC_DIR`, which see.
""" User.music

@doc """
    pictures(parts...) -> String

Join the pictures directory with zero or more path components (`parts`).

The pictures directory is based on the variable `BaseDirs.PICTURES_DIR`, which see.
""" User.pictures

@doc """
    videos(parts...) -> String

Join the videos directory with zero or more path components (`parts`).

The videos directory is based on the variable `BaseDirs.VIDEOS_DIR`, which see.
""" User.videos

@doc """
    templates(parts...) -> String

Join the templates directory with zero or more path components (`parts`).

The templates directory is based on the variable `BaseDirs.TEMPLATES_DIR`, which see.
""" User.templates

@doc """
    public(parts...) -> String

Join the public directory with zero or more path components (`parts`).

The public directory is based on the variable `BaseDirs.PUBLIC_DIR`, which see.
""" User.public

# ---------

@doc """
**BaseDirs.System**

This module contains accessor functions for system directories.

### Base directory accessors

The functions [`data`](@ref BaseDirs.System.data) and [`config`](@ref
BaseDirs.System.config) produce a `Vector{String}`.

### Other accessors

The functions [`fonts`](@ref BaseDirs.System.fonts) and [`applications`](@ref
BaseDirs.System.applications) produce a `Vector{String}`.

See also: [`BaseDirs.User`](@ref).
""" System

@doc """
$(Internals.accessordoc(:data, :DATA_DIRS, name="system configuration"))
""" System.data

@doc """
$(Internals.accessordoc(:config, :CONFIG_DIRS, name="system configuration"))
""" System.config

@doc """
$(Internals.accessordoc(:fonts, name="system fonts", plural=true))
""" System.fonts

@doc """
$(Internals.accessordoc(:applications, name="system applications", plural=true))
""" System.applications

# ---------

@doc """
$(Internals.accessordoc(:data, [:DATA_HOME, :DATA_DIRS], name="user and system configuration"))
""" data

@doc """
$(Internals.accessordoc(:config, [:CONFIG_HOME, :CONFIG_DIRS], name="user and system configuration"))
""" config

@doc """
$(Internals.accessordoc(:fonts, name="user and system fonts", plural=true))
""" fonts

@doc """
$(Internals.accessordoc(:applications, name="user and system applications", plural=true))
""" applications
