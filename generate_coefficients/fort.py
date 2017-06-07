from sympy import Eq, Indexed, IndexedBase, postorder_traversal, fcode, \
    Symbol, Function, Wild, Derivative, MatrixSymbol

def function_derivatives(eq, functions, fval = None):
    functions_extra = []
    for ifun,f in zip(range(len(functions)),functions):
        nargs = len(f.args)
        ff = f.func
        w = [ Wild('w{}'.format(i)) for i in range(nargs) ]
        for i in range(nargs):
            replace_what = Derivative(ff(*w),w[i])
            if fval != None:
                eq = eq.replace( replace_what, fval[i+1,ifun] )
            else:
                f2 = Function('{}_{}'.format(str(ff),i+1))(*w)
                functions_extra.append(f2)
                eq = eq.replace( replace_what, f2 )
        if fval != None:
            eq = eq.replace( ff(*w), fval[0,ifun] )
    return eq.doit(), functions_extra

class FortranProcedure:

    def __init__(self, name, expr_list,
                 extern_functions = {},
                 extern_variables = {},
                 arguments = [],
                 iso_c_binding = False,
                 kind = 8,
                 pure = True):

        self.name = name
        self.real_kind = kind

        assert type(extern_functions) == dict
        assert type(extern_variables) == dict
        assert type(arguments) == list

        self.iso_c_binding = iso_c_binding
        self.pure = pure

        def fixglobs(ex):
            for k,v in extern_variables.items():
                ex = ex.subs(k,Symbol(v))
            return ex

        expr_list_final = [ Eq( fixglobs(ex.lhs), fixglobs(ex.rhs) ) for ex in expr_list ]
        self.expressions = expr_list_final
        self.extern_functions = extern_functions

        def extract_variables(expr):
            # expand the expression tree
            expr_tree = list(postorder_traversal(expr))
            # extract arrays
            arrays = set(filter(lambda x: type(x) == IndexedBase or type(x) == MatrixSymbol, expr_tree))
            # we collect all integer that appear as array dimensions
            array_dimensions = set()
            for a in arrays:
                array_dimensions.update(
                    filter(lambda x: x.is_integer and x.is_Symbol, a.shape)
                )
            # here are symbols corresponding to arrays (arrays appear twice:
            # as IndexedBase and as Symbol)
            array_symbols = set( a.args[0] for a in arrays )
            # extract symbols and remove array symbols
            variables = set(filter(lambda x: x.is_Symbol, expr_tree)) - array_symbols
            return (variables - array_dimensions, arrays, array_dimensions)

        # left and right hand side variables and arrays
        var_lhs, arr_lhs = set(), set()
        var_rhs, arr_rhs = set(), set()
        # any symbols that appear in dimensions (such as 3xN array)
        dims = set()

        # iterate through code lines
        for ex in expr_list_final:
            # extract variables from the LHS
            lhs = extract_variables(ex.lhs)
            var_lhs.update(lhs[0])
            arr_lhs.update(lhs[1])
            dims.update(lhs[2])
            # extract variables from the RHS
            rhs = extract_variables(ex.rhs)
            var_rhs.update(rhs[0])
            arr_rhs.update(rhs[1])
            dims.update(rhs[2])

        # convert the given list of global variables to a set of symbols
        extern_variables_set = set(Symbol(s) for s in extern_variables.values())
        # remove global variables from variables to be declared
        var_lhs -= extern_variables_set
        var_rhs -= extern_variables_set
        arr_lhs -= extern_variables_set
        arr_rhs -= extern_variables_set
        dims    -= extern_variables_set

        # divide arrays and variables basen on whether they occur on the RHS or the LHS
        self.var_in, self.var_inout, self.var_out \
            = var_rhs-var_lhs, var_lhs&var_rhs, var_lhs-var_rhs
        self.arr_in, self.arr_inout, self.arr_out \
            = arr_rhs-arr_lhs, arr_lhs&arr_rhs, arr_lhs-arr_rhs
        self.dims = dims

        self.var_all = self.var_in | self.var_out | self.var_inout
        self.arr_all = self.arr_in | self.arr_out | self.arr_inout

        argument_set = set(arguments)
        for arg in arguments:
            if arg in extern_variables_set:
                arguments.remove(arg)
        self.arguments = arguments + \
                list(self.var_in - argument_set) + \
                list(self.arr_in - argument_set) + \
                list(self.var_inout - argument_set) + \
                list(self.arr_inout - argument_set) + \
                list(self.var_out - argument_set) + \
                list(self.arr_out - argument_set) + \
                list(self.dims - argument_set)
        self.user_arguments = argument_set - (self.var_all | self.arr_all | self.dims) - extern_variables_set


    def __str__(self):

        output = []

        # helper functions to convert set of symbols to a list of
        # uppercase/lowercase strings
        s2s = lambda(S): [ str(s).lower() if s.is_Symbol else str(s).upper() for s in S ]

        output.append('{pure} SUBROUTINE {name}({arglist}) {bindc}'.format(
            name = self.name.upper(),
            arglist = ', '.join(s2s(self.arguments)),
            bindc = 'BIND(C)' if self.iso_c_binding else '',
            pure = 'PURE' if self.pure else '',
        ))

        if self.iso_c_binding: output.append('USE ISO_C_BINDING')

        # this helper prints the fortran declaration of rank-0 variable
        def fort_variables(inset, intent = 'in'):
            l = []
            for v in inset:
                l.append("{vartype}, INTENT({intent}){byvalue} :: {var}".format(
                        var = str(v).lower(),
                        intent = intent,
                        vartype = 'INTEGER' if v.is_integer else 'real({kind})'.format(
                            kind = 'c_double' if self.iso_c_binding else self.real_kind,
                        ),
                        byvalue = ', VALUE' if self.iso_c_binding and intent == 'in' else '',
                ))
            return l

        # is helper prints the declaration of an array
        def fort_matrices(matrices, intent = 'in'):
            l = []
            for m in matrices:
                l.append("real(fp), INTENT({intent}), DIMENSION({dim}) :: {mx}".format(
                    mx = str(m).upper(),
                    intent = intent,
                    dim = ", ".join([ "{}:{}".format(0,d-1) for d in m.shape ])
                ))
            return l

        for d in self.dims:
            output.append("INTEGER, INTENT(in){byvalue} :: {name} ! array dimension".format(
                name = d,
                byvalue = ', VALUE' if self.iso_c_binding else '',
            ))
        output.extend(fort_variables(self.var_in, "in"))
        output.extend(fort_matrices(self.arr_in, "in"))
        output.extend(fort_variables(self.var_inout, "inout"))
        output.extend(fort_matrices(self.arr_inout, "inout"))
        output.extend(fort_variables(self.var_out, "out"))
        output.extend(fort_matrices(self.arr_out, "out"))
        output.extend(fort_variables(self.user_arguments, "in"))

        # iterate through code lines and append the output of fortran code
        # generator from sympy
        for ex in self.expressions:
            output.append(fcode(ex.rhs, assign_to = ex.lhs,
                        contract = False,
                        source_format='free', standard = 2003,
                        user_functions = {str(k): v for k,v in self.extern_functions.items()}))

        # subroutine end
        output.append('END SUBROUTINE {}'.format(self.name.upper()))
        return "\n".join(output)
