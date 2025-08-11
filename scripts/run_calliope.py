import calliope
import xarray as xr
from calliope.backend.helper_functions import ParsingHelperFunction

# Just loading this class in creates the hook within calliope to make this math helper function (sum_next_n) available.
class SumNextN(ParsingHelperFunction):
    """Sum the next N items in an array. Works best for ordered arrays (datetime, integer)."""

    #:
    NAME = "sum_next_n"
    #:
    ALLOWED_IN = ["expression"]

    def as_math_string(self, array: str, *, over: str, N: int) -> str:  # noqa: D102, override
        overstring = self._instr(over)
        # FIXME: add N
        return rf"\sum\limits_{{{overstring}}} ({array})"

    def as_array(self, array: xr.DataArray, over: str, N: int) -> xr.DataArray:
        """Sum values from current up to N from current on the dimension `over`.

        Args:
            array (xr.DataArray): Math component array.
            over (str): Dimension over which to sum
            N (int): number of items beyond the current value to sum from

        Returns:
            xr.DataArray:
                Returns the input array with the condition applied,
                including having been broadcast across any new dimensions provided by the condition.

        Examples:
            One common use-case is to collate N timesteps beyond a given timestep to apply a constraint to it(e.g., demand must be less than X in the next 24 hours)
        """
        results = []
        for i in range(len(self._input_data.coords[over])):
            results.append(
                array.isel(**{over: slice(i, i + int(N))}).sum(over, min_count=1)
            )
        final_array = xr.concat(
            results, dim=self._input_data.coords[over]
        ).broadcast_like(array)
        return final_array


def create_dsr_constraint(response_hrs: int, resolution_hrs: int) -> dict:
    """Create demand-side response math that allows shifting demands within a specified rolling time range.

    Args:
        response_hrs (int): Number of rolling hours over which to allow DSR
        resolution_hrs (int): Model time resolution in hours

    Raises:
        ValueError: response_hrs / resolution_hrs must be an integer

    Returns:
        dict: Calliope math constraint definition
    """
    if response_hrs % resolution_hrs != 0:
        raise ValueError(
            f"Cannot run a model with {response_hrs} hr demand-side response when running a model with a {resolution_hrs} hr resolution of"
        )
    n_timesteps = response_hrs // resolution_hrs
    dsr = {
        f"demand_side_response_{n_timesteps}hr": {
            "description": """
            Demand side response, allowing X% demand to be shifted within a rolling window of 4 hours.
            Will work with any technology that has an inflow and a `sink_use_dsr` timeseries parameter set.
            """,
            "foreach": ["nodes", "techs", "carriers", "timesteps"],
            # Do not set the constraint in the final N timesteps otherwise it over-constrains the model
            # such that `sink_use_dsr` operates the same way as `sink_use_equals` would.
            "where": f"carrier_in AND sink_use_dsr AND timesteps<=get_val_at_index(timesteps=-{n_timesteps:d})",
            "equations": [
                {
                    "expression": f"sum_next_n(flow_in, timesteps, {n_timesteps:d}) == sum_next_n(sink_use_dsr, timesteps, {n_timesteps:d})"
                }
            ],
        }
    }
    return dsr


math = {
    "constraints": create_dsr_constraint(
        int(snakemake.params.response_hrs), int(snakemake.params.resolution_hrs)
    )
}


model = calliope.Model(
    snakemake.input[0], time_resample=f"{snakemake.params.resolution_hrs}h"
)

calliope.set_log_verbosity("INFO")
model.build(add_math_dict=math)
#model.build()
model.solve()

model.to_netcdf(snakemake.output[0])
