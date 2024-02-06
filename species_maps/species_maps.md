# Species map TOML file format

## Example

```toml
[meta]
name = "Map name"
version = 1

[units]
source = "kg"
target = "ug"

[aliases_source]
VAR1 = "asdf1"
VAR2 = "asdf2"
VAR3 = "asdf3"

[aliases_target]
TARGET1 = "my_var_bin1"
TARGET2 = "my_var_bin2"

[species_map]
TARGET1 = { VAR1 = 1.0, VAR2 = 0.5 }
TARGET2 = { VAR2 = 0.5, VAR3 = 1.0 }
```

## `meta` section

| --- | --- | --- |
| Key | Type | Description |
| --- | --- | --- |
| name | string | Name of the species map |
| version | int | Version of the application, set to `1` |

## `units` section

| --- | --- | --- |
| Key | Type | Description |
| --- | --- | --- |
| source | string | Units of the source fields |
| target | string | Units of the target fields |

Possible values are: `kg`, `g`, `mg`, `ug`. If you set both `source` and target to `kg`, then no unit conversion takes place (no-op).

## `aliases_source` and `aliases_target` sections

These sections are optional. They are used to map the names of the source and target fields to the actual names of the fields in the global model and WRF-CHEM, respectively. For example, if the source model has a useful field with a weird name like `ASDF1`, you can map it to a more meaningful name like `DUST1` in the `aliases_source` section. You can use these aliases in the `species_map` section.

## `species_map` section

This section contains the species mapping. Each key is the variable name of WRF-CHEM (or a target alias) and the value should be a dictionary of source variables and their coefficients. The coefficients are used to create the target variable through a linear combination of the source variables. For example, if you want to create `DUST1` in WRF-CHEM through a linear combination of `ASDF1` and `ASDF2` from the global model, you can use the following:

```toml
[species_map]
DUST1 = { ASDF1 = 1.0, ASDF2 = 0.5 }
```

You can use any variable from the source model, even if it's not in the `aliases_source` section. If the variable is not in the `aliases_source` section, the variable name is used as is.
