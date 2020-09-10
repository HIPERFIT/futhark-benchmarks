#!/usr/bin/env python3

import sys
import re

out_file, in_file = sys.argv[1:]

with open(in_file) as f:
    contents = f.read()

start = contents.find('futhark_entry_text_content')
end = contents.find(')', start)
types = re.findall('([^ ]+) \*out\d+,', contents[start:end])
out_vars = ['out{}'.format(i) for i in range(len(types))]

with open(out_file, 'w') as f:
    print('#include <stdio.h>', file=f)
    print('#include "lib/github.com/diku-dk/lys/liblys.h"', file=f)
    print('', file=f)
    if len(types) == 0:
        print('#define UNUSED(x) (void)(x)', file=f)
    print('void build_text(const struct lys_context *ctx, char* dest, size_t dest_len, const char* format, float render_milliseconds, char* **sum_names) {', file=f)
    if len(types) > 0:
        for v, t in zip(out_vars, types):
            print('  union {{ {} val; char* sum_name; }} {};'.format(t, v), file=f)
        print('  FUT_CHECK(ctx->fut, futhark_entry_text_content(ctx->fut, {}, render_milliseconds, ctx->state));'.format(', '.join('&{}.val'.format(v) for v in out_vars)), file=f)
        for v, i in zip(out_vars, range(len(out_vars))):
            print('  if (sum_names[{}] != NULL) {{'.format(i), file=f)
            print('    {v}.sum_name = sum_names[{i}][(int32_t) {v}.val];'.format(v=v, i=i), file=f)
            print('  }', file=f)
        print('  snprintf(dest, dest_len, format, {});'.format(', '.join((s + ('.sum_name' if t == 'int32_t' else '.val')) for s, t in zip(out_vars, types))), file=f)
    else:
        for x in ['ctx', 'render_milliseconds', 'sum_names']:
            print('UNUSED({});'.format(x), file=f)
        print('  snprintf(dest, dest_len, "%s", format);', file=f)
    print('}', file=f)
    print('', file=f)
    print('size_t n_printf_arguments() {', file=f)
    print('  return {};'.format(len(types)), file=f)
    print('}', file=f)
