import numpy as np


def object_einsum(string, *arrays):
    """Simplified object einsum, not as much error checking

    does not support "..." or list input and will see "...", etc. as three times
    an axes identifier, tries normal einsum first!

    NOTE: This is untested, and not fast, but object type is
    never really fast anyway...
    """
    try:
        return np.einsum(string, *arrays)
    except TypeError:
        pass

    s = string.split('->')
    in_op = s[0].split(',')
    out_op = None if len(s) == 1 else s[1].replace(' ', '')

    in_op = [axes.replace(' ', '') for axes in in_op]
    all_axes = set()
    repeated_axes = set()

    for axes in in_op:
        list(repeated_axes.update(ax) for ax in axes if ax in all_axes)
        all_axes.update(axes)

    if out_op is None:
        out_op = set(sorted(all_axes))
        list(out_op.discard(rep_ax) for rep_ax in repeated_axes)
    else:
        all_axes.update(out_op)

    perm_dict = {_[1]: _[0] for _ in enumerate(all_axes)}

    dims = len(perm_dict)
    op_axes = []
    for axes in (in_op + list((out_op,))):
        op = [-1] * dims
        for i, ax in enumerate(axes):
            op[perm_dict[ax]] = i
        op_axes.append(op)

    op_flags = [('readonly',)] * len(in_op) + [('readwrite', 'allocate')]
    dtypes = [np.object_] * (len(in_op) + 1)  # cast all to object

    nditer = np.nditer(arrays + (None,), op_axes=op_axes,
                       flags=['buffered', 'delay_bufalloc', 'reduce_ok', 'grow_inner', 'refs_ok'], op_dtypes=dtypes,
                       op_flags=op_flags)

    nditer.operands[-1][...] = 0
    nditer.reset()

    for vals in nditer:
        out = vals[-1]
        prod = np.copy(vals[0])
        for value in vals[1:-1]:
            prod *= value
        out += prod

    return nditer.operands[-1]