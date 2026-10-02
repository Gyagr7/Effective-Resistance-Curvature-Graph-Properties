/*
 * rn.c -- recognise RN / RP / SRN graphs.
 *
 * AI disclosure.
 * Claude Opus 5 (Anthropic) wrote this file in its entirety: the formulation of the program, the structural shortcuts, the separation routine, the exact rational certification, the I/O, and these comments.
 * Hailey Jay Garcia directed and supervised that work, reviewed the result, and is responsible for the correctness of the mathematics it rests on and of the answers it reports.
 *
 * Citation keys are those of ../references.bib.
 *
 * A graph G is resistance nonnegative (RN) if some positive conductance function c makes the resistance curvature p_v(c) = 1 - (1/2) sum_{u~v} r_uv(c) nonnegative at every vertex.
 * It is resistance positive (RP) if some c makes that curvature strictly positive everywhere, and strictly resistance nonnegative (SRN) if it is RN but not RP.
 * Devriendt's polytope characterisation (devriendt2025) turns this into
 *
 *     G is RN  <=>  P(G)^o INTERSECT { x : d_v(x) <= 2 for all v }  !=  empty
 *     G is RP  <=>  P(G)^o INTERSECT { x : d_v(x) <  2 for all v }  !=  empty
 *
 * where P(G) is the spanning tree polytope, P(G)^o its relative interior, and d_v(x) = sum_{e ni v} x_e.
 * (guo2026lp, Theorem 2.1) makes both decidable by linear programming.
 *
 * The delicate point is the relative interior, which one does not get by making the inequalities strict the moment G has a bridge.
 * (guo2026lp) handles it by measuring slack against an explicit interior point, the uniform spanning tree marginal vector.
 * We take the other route and remove the difficulty structurally, by a lemma that leaves the program a single case to decide.
 *
 *   Lemma.  A connected graph that is not 2-connected is RN exactly when it is a path.
 *
 *   Proof.  Let v be a cut vertex of G, lying in the blocks B_1, ..., B_k, k >= 2.
 *   Spanning trees of G restrict to spanning trees of its blocks and conversely, so P(G) is the product of the P(B_i) and d_v(x) is the sum of the d_v^{B_i}(x).
 *   Every spanning tree of B_i uses an edge at v, so d_v^{B_i}(x) >= 1 on P(B_i), with equality for all x only when B_i is a bridge.
 *   When B_i is 2-connected, v has two neighbours in it and some spanning tree of B_i uses both of those edges, so d_v^{B_i}(x) > 1 on P(B_i)^o.
 *   Hence d_v(x) > 2 on P(G)^o unless k = 2 and both blocks are bridges.
 *   So in an RN graph every block at a cut vertex is a bridge and every cut vertex lies in exactly two of them: G is a tree of maximum degree 2, a path.
 *   Conversely P_n is RN at x = 1, and RP exactly when n <= 2.  []
 *
 * A path is therefore answered outright, and so is every other graph with a cut vertex, which leaves the program only 2-connected graphs.
 * A counting bound over vertex cuts settles more of them for nothing.
 *
 *   Lemma.  Let S be a nonempty proper set of vertices and c the number of components of G - S.  Then
 *       max_{v in S} d_v(x) >= 1 + (c - 1)/|S|   for every x in P(G).
 *   So c >= |S| + 1 makes G not RP, and c >= |S| + 2 makes G not RN.
 *
 *   Proof.  Write C_1, ..., C_c for the components of G - S.  Every edge lies inside S, inside one C_i, or across, so
 *       sum_{v in S} d_v(x) = 2 x(E[S]) + x(dS) = x(E[S]) + (n - 1 - sum_i x(E[C_i]))
 *   using x(E) = n - 1.  The rank inequalities give x(E[C_i]) <= |C_i| - 1 and x >= 0 gives x(E[S]) >= 0, so
 *       sum_{v in S} d_v(x) >= (n - 1) - (n - |S| - c) = |S| + c - 1,
 *   and the claim follows by averaging over S.  []
 *
 * Both halves are known, and the certifier only needs the cut that witnesses them: the RP half is the 1-toughness of resistance positive graphs (devriendt2025, by the toughness argument of chvatal1973), and the RN half is the sharpening in (garcia2026tough).
 * Taking S to be a vertex cover, where every C_i is a single vertex, gives the bipartite corner: either part of a bipartite graph has c = n - |S|, so parts differing by one rule out RP and parts differing by two or more rule out RN.
 * The worst cut is NP-hard to find, but any cut is a certificate on its own, so the search below is a fixed, deterministic list of candidates: it is sound wherever it fires and costs one sweep where it does not.
 * A cut with c >= |S| + 2 therefore answers the graph outright, and one with c = |S| + 1 skips the RP program and runs only the RN one, which is where the saving is: the graphs that are RN but not RP are exactly the ones that would otherwise pay for both programs.
 * There P(G) carries the one equality x(E) = n - 1 and no implicit one, so there the relative interior really is the strict system.
 * That leaves one program with integer data, and integer data is what makes exact rational arithmetic practical:
 *
 *     max t   s.t.   x(E) = n - 1
 *                    x_e >= t                                   (every edge e)
 *                    x(E[S]) + (|S|-1) t <= |S|-1     (S a proper nonempty set
 *                                                            of vertices of G)
 *                    d_v(x) + s t <= 2                              (every v)
 *
 * with s = 0 for RN and s = 1 for RP, x >= 0 and -n <= t <= 1.
 * The program is always feasible, and the property holds exactly when the optimum is positive.
 * An optimal x is then a witness outright.
 * A negative RP optimum says more: no point of P(G) has all degrees at most 2, so G is not RN either, and the RN program can be skipped.
 * The (|S|-1) coefficient in place of a flat t keeps singletons from registering as violated.
 * The separation below therefore runs in O(n) flows rather than O(n^3).
 *
 * The rank inequalities number 2^n, so they are generated on demand.
 * Since
 *
 *     (1-t)|S| - x(E[S]) = sum_{v in S} ((1-t) - d_v(x)/2) + x(dS)/2
 *
 * is a cut function, a most violated S is a minimum s-t cut, and one round of separation is 2(n-1) maximum flows.
 * Fix a vertex v0 and, for each other u, force u out of S and then into it; between them those cover every proper nonempty S.
 * That cost is what puts graphs on a hundred vertices in reach.
 *
 * With --exact the answer is then certified in rational arithmetic.
 * The tight rows at the floating point optimum are selected down to a square nonsingular system, which is solved exactly both ways.
 * The forward solve gives the vertex, checked against every constraint including an exact rerun of the separation.
 * The transposed solve gives the dual, checked for sign feasibility, and it yields an exact upper bound on the optimum.
 * A positive certified vertex proves the property, and a nonpositive certified dual bound disproves it.
 * Neither is assumed: when a certificate does not come out, the output says so.
 *
 * Input is JSON or graph6, and output is JSON.
 * See README.md.
 *
 * Build:  make        (needs GLPK and GMP: -lglpk -lgmp)
 */

#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <glpk.h>
#include <gmp.h>

/* Sign of an LP optimum.
   Both optima are ratios of subdeterminants of an integer matrix of order |E|+1, so a true positive is never this small.
   Under --exact the sign is settled in rational arithmetic, and this tolerance is only a hint about where to look. */
#define OPT_TOL  1e-9
/* Slack allowed of a row when it is checked after the fact, and when it is read as tight for rebuilding the vertex exactly. */
#define FEAS_TOL 1e-7
/* A rank inequality counts as violated only beyond this. */
#define CUT_TOL  1e-7
#define MAX_ROUNDS 2000

static void die(const char *msg)
{
	fprintf(stderr, "rn: %s\n", msg);
	exit(2);
}

static void *xmalloc(size_t k)
{
	void *p = malloc(k ? k : 1);
	if (!p)
		die("out of memory");
	return p;
}

static void *xcalloc(size_t k, size_t s)
{
	void *p = calloc(k ? k : 1, s ? s : 1);
	if (!p)
		die("out of memory");
	return p;
}

static void *xrealloc(void *p, size_t k)
{
	void *q = realloc(p, k ? k : 1);
	if (!q)
		die("out of memory");
	return q;
}

static char *xstrdup(const char *s)
{
	size_t k = strlen(s) + 1;
	char *t = xmalloc(k);
	memcpy(t, s, k);
	return t;
}

/* ------------------------------------------------------------------ JSON in */

/* A lax reader: commas and colons count as whitespace, so it accepts more than JSON proper.
   That is deliberate, since this reads input we are handed rather than validating it. */

typedef struct JV JV;
struct JV {
	enum { JNULL, JBOOL, JNUM, JSTR, JARR, JOBJ } t;
	double num;		/* JNUM, and 0/1 for JBOOL */
	char *str;		/* JSTR */
	JV **kid;		/* JARR, JOBJ */
	char **key;		/* JOBJ */
	int nk, cap;
};

typedef struct {
	const char *p;
	const char *err;
} JP;

static JV *jvalue(JP *j);

static JV *jnew(int t)
{
	JV *v = xcalloc(1, sizeof *v);
	v->t = t;
	return v;
}

static void jfree(JV *v)
{
	int i;
	if (!v)
		return;
	for (i = 0; i < v->nk; i++) {
		jfree(v->kid[i]);
		if (v->key)
			free(v->key[i]);
	}
	free(v->kid);
	free(v->key);
	free(v->str);
	free(v);
}

static void jpush(JV *v, char *key, JV *kid)
{
	if (v->nk == v->cap) {
		v->cap = v->cap ? 2 * v->cap : 8;
		v->kid = xrealloc(v->kid, (size_t)v->cap * sizeof *v->kid);
		if (v->t == JOBJ)
			v->key = xrealloc(v->key, (size_t)v->cap * sizeof *v->key);
	}
	if (v->t == JOBJ)
		v->key[v->nk] = key;
	v->kid[v->nk++] = kid;
}

static void jskip(JP *j)
{
	while (*j->p && (isspace((unsigned char)*j->p) || *j->p == ',' || *j->p == ':'))
		j->p++;
}

/* Append one code point as UTF-8. */
static void utf8(char **out, size_t *len, size_t *cap, unsigned cp)
{
	char buf[4];
	int k = 0;
	if (cp < 0x80) {
		buf[k++] = (char)cp;
	} else if (cp < 0x800) {
		buf[k++] = (char)(0xC0 | (cp >> 6));
		buf[k++] = (char)(0x80 | (cp & 0x3F));
	} else {
		buf[k++] = (char)(0xE0 | (cp >> 12));
		buf[k++] = (char)(0x80 | ((cp >> 6) & 0x3F));
		buf[k++] = (char)(0x80 | (cp & 0x3F));
	}
	if (*len + (size_t)k + 1 > *cap) {
		*cap = 2 * (*cap) + (size_t)k + 16;
		*out = xrealloc(*out, *cap);
	}
	memcpy(*out + *len, buf, (size_t)k);
	*len += (size_t)k;
}

static char *jstring(JP *j)
{
	char *out = NULL;
	size_t len = 0, cap = 0;
	if (*j->p != '"') {
		j->err = "expected a string";
		return NULL;
	}
	j->p++;
	out = xmalloc(cap = 16);
	while (*j->p && *j->p != '"') {
		unsigned cp;
		if (*j->p != '\\') {
			utf8(&out, &len, &cap, (unsigned char)*j->p++);
			continue;
		}
		j->p++;
		switch (*j->p) {
		case 'n': cp = '\n'; j->p++; break;
		case 't': cp = '\t'; j->p++; break;
		case 'r': cp = '\r'; j->p++; break;
		case 'b': cp = '\b'; j->p++; break;
		case 'f': cp = '\f'; j->p++; break;
		case 'u': {
			char hex[5] = { 0 };
			int i;
			j->p++;
			for (i = 0; i < 4 && j->p[i]; i++)
				hex[i] = j->p[i];
			cp = (unsigned)strtoul(hex, NULL, 16);
			j->p += i;
			break;
		}
		default: cp = (unsigned char)*j->p; if (*j->p) j->p++; break;
		}
		utf8(&out, &len, &cap, cp);
	}
	if (*j->p != '"') {
		j->err = "unterminated string";
		free(out);
		return NULL;
	}
	j->p++;
	out[len] = 0;
	return out;
}

static JV *jvalue(JP *j)
{
	JV *v;
	jskip(j);
	switch (*j->p) {
	case 0:
		j->err = "unexpected end of input";
		return NULL;
	case '{':
		j->p++;
		v = jnew(JOBJ);
		for (;;) {
			char *k;
			JV *c;
			jskip(j);
			if (*j->p == '}') { j->p++; return v; }
			if (!*j->p) { j->err = "unterminated object"; jfree(v); return NULL; }
			if (!(k = jstring(j))) { jfree(v); return NULL; }
			if (!(c = jvalue(j))) { free(k); jfree(v); return NULL; }
			jpush(v, k, c);
		}
	case '[':
		j->p++;
		v = jnew(JARR);
		for (;;) {
			JV *c;
			jskip(j);
			if (*j->p == ']') { j->p++; return v; }
			if (!*j->p) { j->err = "unterminated array"; jfree(v); return NULL; }
			if (!(c = jvalue(j))) { jfree(v); return NULL; }
			jpush(v, NULL, c);
		}
	case '"':
		v = jnew(JSTR);
		if (!(v->str = jstring(j))) { jfree(v); return NULL; }
		return v;
	case 't':
		j->p += 4; v = jnew(JBOOL); v->num = 1; return v;
	case 'f':
		j->p += 5; v = jnew(JBOOL); v->num = 0; return v;
	case 'n':
		j->p += 4; return jnew(JNULL);
	default: {
		char *end;
		double d = strtod(j->p, &end);
		if (end == j->p) { j->err = "expected a value"; return NULL; }
		j->p = end;
		v = jnew(JNUM);
		v->num = d;
		return v;
	}
	}
}

static JV *jget(const JV *o, const char *k)
{
	int i;
	if (!o || o->t != JOBJ)
		return NULL;
	for (i = 0; i < o->nk; i++)
		if (!strcmp(o->key[i], k))
			return o->kid[i];
	return NULL;
}

/* ----------------------------------------------------------------- JSON out */

static void jsquote(FILE *f, const char *s)
{
	fputc('"', f);
	for (; *s; s++) {
		unsigned char c = (unsigned char)*s;
		switch (c) {
		case '"':  fputs("\\\"", f); break;
		case '\\': fputs("\\\\", f); break;
		case '\n': fputs("\\n", f); break;
		case '\r': fputs("\\r", f); break;
		case '\t': fputs("\\t", f); break;
		default:
			if (c < 0x20)
				fprintf(f, "\\u%04x", c);
			else
				fputc(c, f);
		}
	}
	fputc('"', f);
}

/* Numbers print short but round-trippable. */
static void jsnum(FILE *f, double x)
{
	if (x == 0.0)
		x = 0.0;		/* never print -0 */
	if (x == floor(x) && fabs(x) < 1e15)
		fprintf(f, "%.0f", x);
	else
		fprintf(f, "%.17g", x);
}

/* -------------------------------------------------------------------- graph */

typedef struct {
	char *name;
	int n, m;
	char **vname;		/* NULL when the vertices were unlabelled */
	int *eu, *ev;		/* endpoints, 0-based, eu < ev */
	char *note;		/* set when the input had to be repaired */
} Graph;

static void gfree(Graph *g)
{
	int i;
	if (!g)
		return;
	if (g->vname) {
		for (i = 0; i < g->n; i++)
			free(g->vname[i]);
		free(g->vname);
	}
	free(g->name);
	free(g->note);
	free(g->eu);
	free(g->ev);
}

static int edge_cmp(const void *a, const void *b)
{
	const int *x = a, *y = b;
	if (x[0] != y[0])
		return x[0] < y[0] ? -1 : 1;
	return x[1] == y[1] ? 0 : (x[1] < y[1] ? -1 : 1);
}

/* Sort the edges and drop loops and repeats.
   Return the number of edges dropped. */
static int normalise(Graph *g)
{
	int *buf, i, k = 0, dropped = 0;
	if (g->m == 0)
		return 0;
	buf = xmalloc((size_t)g->m * 2 * sizeof *buf);
	for (i = 0; i < g->m; i++) {
		int u = g->eu[i], v = g->ev[i];
		buf[2 * i] = u < v ? u : v;
		buf[2 * i + 1] = u < v ? v : u;
	}
	qsort(buf, (size_t)g->m, 2 * sizeof *buf, edge_cmp);
	for (i = 0; i < g->m; i++) {
		if (buf[2 * i] == buf[2 * i + 1])
			continue;	/* loop */
		if (k && buf[2 * k - 2] == buf[2 * i] && buf[2 * k - 1] == buf[2 * i + 1])
			continue;	/* repeat */
		buf[2 * k] = buf[2 * i];
		buf[2 * k + 1] = buf[2 * i + 1];
		k++;
	}
	dropped = g->m - k;
	for (i = 0; i < k; i++) {
		g->eu[i] = buf[2 * i];
		g->ev[i] = buf[2 * i + 1];
	}
	g->m = k;
	free(buf);
	return dropped;
}

static int connected(const Graph *g)
{
	int *stack, *seen, top = 0, cnt = 1, i;
	int *head, *next, *to;
	if (g->n <= 1)
		return 1;
	head = xmalloc((size_t)g->n * sizeof *head);
	next = xmalloc((size_t)g->m * 2 * sizeof *next);
	to = xmalloc((size_t)g->m * 2 * sizeof *to);
	for (i = 0; i < g->n; i++)
		head[i] = -1;
	for (i = 0; i < g->m; i++) {
		to[2 * i] = g->ev[i]; next[2 * i] = head[g->eu[i]]; head[g->eu[i]] = 2 * i;
		to[2 * i + 1] = g->eu[i]; next[2 * i + 1] = head[g->ev[i]]; head[g->ev[i]] = 2 * i + 1;
	}
	seen = xcalloc((size_t)g->n, sizeof *seen);
	stack = xmalloc((size_t)g->n * sizeof *stack);
	stack[top++] = 0;
	seen[0] = 1;
	while (top) {
		int v = stack[--top], a;
		for (a = head[v]; a >= 0; a = next[a])
			if (!seen[to[a]]) {
				seen[to[a]] = 1;
				cnt++;
				stack[top++] = to[a];
			}
	}
	free(head); free(next); free(to); free(seen); free(stack);
	return cnt == g->n;
}

/* ------------------------------------------------------------------- graph6 */

/* Read the order from graph6's N(n): one byte under 63, otherwise `~' and three more, or `~~' and six.
   Advance *p past it, and return -1 on a bad header. */
static long g6_order(const unsigned char **p)
{
	const unsigned char *q = *p;
	int k, i;
	long n = 0;
	if (!*q || *q < 63 || *q > 126)
		return -1;
	if (*q != 126) { *p = q + 1; return *q - 63; }
	q++;
	k = (*q == 126) ? (q++, 6) : 3;
	for (i = 0; i < k; i++) {
		if (q[i] < 63 || q[i] > 126)
			return -1;
		n = (n << 6) | (q[i] - 63);
	}
	*p = q + k;
	return n;
}

static int g6_decode(const char *s, Graph *g)
{
	int i, j, k, bit = 0;
	long n;
	const unsigned char *p = (const unsigned char *)s;
	size_t need;
	n = g6_order(&p);
	if (n < 0 || n > 100000)
		return 0;
	need = ((size_t)n * (size_t)(n - 1) / 2 + 5) / 6;
	if (strlen((const char *)p) < need)
		return 0;
	for (i = 0; i < (int)need; i++)
		if (p[i] < 63 || p[i] > 126)
			return 0;
	g->n = (int)n;
	g->m = 0;
	g->eu = xmalloc((size_t)n * (size_t)n * sizeof *g->eu / 2 + sizeof *g->eu);
	g->ev = xmalloc((size_t)n * (size_t)n * sizeof *g->ev / 2 + sizeof *g->ev);
	for (j = 1; j < n; j++)
		for (i = 0; i < j; i++) {
			k = bit++;
			if ((p[k / 6] - 63) >> (5 - k % 6) & 1) {
				g->eu[g->m] = i;
				g->ev[g->m] = j;
				g->m++;
			}
		}
	return 1;
}

/* The caller frees the result. */
static char *g6_encode(const Graph *g)
{
	int n = g->n, i, j, k, hd;
	size_t bits = (size_t)n * (size_t)(n - 1) / 2, nb = (bits + 5) / 6;
	unsigned char *adj;
	char *out;
	if (n > 258047)
		return NULL;
	adj = xcalloc((size_t)n * (size_t)n, 1);
	for (i = 0; i < g->m; i++)
		adj[(size_t)g->eu[i] * n + g->ev[i]] = adj[(size_t)g->ev[i] * n + g->eu[i]] = 1;
	hd = n < 63 ? 1 : 4;
	out = xcalloc(nb + (size_t)hd + 1, 1);
	if (hd == 1) {
		out[0] = (char)(n + 63);
	} else {
		out[0] = 126;
		for (i = 0; i < 3; i++)
			out[1 + i] = (char)(((n >> (6 * (2 - i))) & 63) + 63);
	}
	for (k = 0; k < (int)nb; k++)
		out[hd + k] = 63;
	k = 0;
	for (j = 1; j < n; j++)
		for (i = 0; i < j; i++) {
			if (adj[(size_t)i * n + j])
				out[hd + k / 6] = (char)(out[hd + k / 6] + (1 << (5 - k % 6)));
			k++;
		}
	free(adj);
	return out;
}

/* ------------------------------------------------- building a graph from JSON */

/* Vertex names are interned in a plain array, since n is small. */
typedef struct {
	char **name;
	int n, cap;
} Names;

static int intern(Names *N, const char *s)
{
	int i;
	for (i = 0; i < N->n; i++)
		if (!strcmp(N->name[i], s))
			return i;
	if (N->n == N->cap) {
		N->cap = N->cap ? 2 * N->cap : 16;
		N->name = xrealloc(N->name, (size_t)N->cap * sizeof *N->name);
	}
	N->name[N->n] = xstrdup(s);
	return N->n++;
}

/* A vertex reference is a string or a number. */
static int vref(Names *N, const JV *v, int *numeric)
{
	char buf[64];
	if (v->t == JSTR)
		{ *numeric = 0; return intern(N, v->str); }
	if (v->t == JNUM) {
		snprintf(buf, sizeof buf, "%.0f", v->num);
		return intern(N, buf);
	}
	*numeric = -1;
	return 0;
}

/* Return 0 and set *err on a malformed entry. */
static int graph_from_json(const JV *o, Graph *g, const char **err)
{
	Names N = { 0 };
	const JV *jv, *je;
	int numeric = 1, i, dropped, want_n = -1;

	memset(g, 0, sizeof *g);
	*err = NULL;

	if (o->t == JSTR) {		/* a bare graph6 string */
		if (!g6_decode(o->str, g)) { *err = "bad graph6 string"; return 0; }
		return 1;
	}
	if (o->t != JOBJ) { *err = "expected an object"; return 0; }

	if ((jv = jget(o, "name")) && jv->t == JSTR)
		g->name = xstrdup(jv->str);

	if ((jv = jget(o, "graph6")) && jv->t == JSTR) {
		if (!g6_decode(jv->str, g)) { *err = "bad graph6 string"; return 0; }
		return 1;
	}

	if ((jv = jget(o, "n")) && jv->t == JNUM)
		want_n = (int)jv->num;

	/* Declared vertices first, so that their order is respected. */
	if ((jv = jget(o, "vertices")) && jv->t == JARR)
		for (i = 0; i < jv->nk; i++)
			vref(&N, jv->kid[i], &numeric);
	else if (want_n > 0)
		for (i = 0; i < want_n; i++) {
			char buf[32];
			snprintf(buf, sizeof buf, "%d", i);
			intern(&N, buf);
		}

	je = jget(o, "edges");
	if (!je)
		je = jget(o, "E");
	if (!je || je->t != JARR) { *err = "no \"edges\" array"; return 0; }

	g->eu = xmalloc((size_t)(je->nk + 1) * sizeof *g->eu);
	g->ev = xmalloc((size_t)(je->nk + 1) * sizeof *g->ev);
	for (i = 0; i < je->nk; i++) {
		const JV *e = je->kid[i];
		if (e->t != JARR || e->nk != 2) { *err = "an edge is not a pair"; return 0; }
		g->eu[i] = vref(&N, e->kid[0], &numeric);
		g->ev[i] = vref(&N, e->kid[1], &numeric);
		if (numeric < 0) { *err = "an endpoint is neither a string nor a number"; return 0; }
	}
	g->m = je->nk;
	g->n = N.n;
	if (want_n > g->n)
		g->n = want_n;

	if (numeric) {		/* every label was a number: keep them implicit */
		for (i = 0; i < N.n; i++)
			free(N.name[i]);
		free(N.name);
		g->vname = NULL;
	} else {
		g->vname = xmalloc((size_t)g->n * sizeof *g->vname);
		for (i = 0; i < N.n; i++)
			g->vname[i] = N.name[i];
		for (i = N.n; i < g->n; i++) {	/* padded by "n" */
			char buf[32];
			snprintf(buf, sizeof buf, "%d", i);
			g->vname[i] = xstrdup(buf);
		}
		free(N.name);
	}

	dropped = normalise(g);
	if (dropped) {
		char buf[96];
		snprintf(buf, sizeof buf, "%d loop or repeated edge%s dropped",
			 dropped, dropped == 1 ? "" : "s");
		g->note = xstrdup(buf);
	}
	return 1;
}

/* ------------------------------------------------- the structural shortcut */

/*
 * Decide 2-connectivity, and name a cut vertex when there is one.
 * Return 1 when G is 2-connected, and 0 otherwise, setting *cut to a cut vertex, or to -1 when there is none (n <= 2).
 * G is already known to be connected, so one scan from vertex 0 sees everything.
 * This is Tarjan's low-link scan, iterative, and it keeps no edge stack: by the lemma at the head of the file the blocks themselves are never needed.
 */
static int two_connected(const Graph *g, int *cut)
{
	int n = g->n, m = g->m, i, a, timer = 0, top = 0, kids = 0, found = -1;
	int *head, *nxt, *to, *disc, *low, *pe, *it, *vstk;

	*cut = -1;
	if (n < 3)
		return 0;

	head = xmalloc((size_t)n * sizeof *head);
	nxt = xmalloc((size_t)m * 2 * sizeof *nxt);
	to = xmalloc((size_t)m * 2 * sizeof *to);
	for (i = 0; i < n; i++)
		head[i] = -1;
	for (i = 0; i < m; i++) {
		to[2 * i] = g->ev[i]; nxt[2 * i] = head[g->eu[i]]; head[g->eu[i]] = 2 * i;
		to[2 * i + 1] = g->eu[i]; nxt[2 * i + 1] = head[g->ev[i]]; head[g->ev[i]] = 2 * i + 1;
	}
	disc = xcalloc((size_t)n, sizeof *disc);
	low = xmalloc((size_t)n * sizeof *low);
	pe = xmalloc((size_t)n * sizeof *pe);
	it = xmalloc((size_t)n * sizeof *it);
	vstk = xmalloc((size_t)n * sizeof *vstk);

	disc[0] = low[0] = ++timer;
	pe[0] = -1;
	it[0] = head[0];
	vstk[top++] = 0;
	while (top) {
		int v = vstk[top - 1];
		if (it[v] >= 0) {
			int w;
			a = it[v];
			it[v] = nxt[a];
			if (pe[v] >= 0 && (a >> 1) == (pe[v] >> 1))
				continue;		/* the edge we came in on */
			w = to[a];
			if (!disc[w]) {
				if (v == 0)
					kids++;
				disc[w] = low[w] = ++timer;
				pe[w] = a;
				it[w] = head[w];
				vstk[top++] = w;
			} else if (disc[w] < low[v]) {
				low[v] = disc[w];
			}
		} else {
			top--;
			if (top) {
				int u = vstk[top - 1];
				if (low[v] < low[u])
					low[u] = low[v];
				if (u != 0 && low[v] >= disc[u] && found < 0)
					found = u;	/* u cuts v off */
			}
		}
	}
	if (found < 0 && kids > 1)
		found = 0;			/* the root, with two subtrees */

	free(head); free(nxt); free(to);
	free(disc); free(low); free(pe); free(it); free(vstk);
	*cut = found;
	return found < 0;
}

/*
 * Two-colour the connected graph G.
 * Return 1 and set *small to the size of the smaller part, or 0 when G has an odd cycle.
 * When out is not NULL it receives the colouring, since both parts are candidate cuts for the second lemma at the head of the file.
 */
static int bipartition(const Graph *g, int *small, unsigned char *out)
{
	int n = g->n, m = g->m, i, top = 0, black = 0, ok = 1;
	int *head, *nxt, *to, *col, *stk;

	*small = 0;
	head = xmalloc((size_t)n * sizeof *head);
	nxt = xmalloc((size_t)m * 2 * sizeof *nxt);
	to = xmalloc((size_t)m * 2 * sizeof *to);
	for (i = 0; i < n; i++)
		head[i] = -1;
	for (i = 0; i < m; i++) {
		to[2 * i] = g->ev[i]; nxt[2 * i] = head[g->eu[i]]; head[g->eu[i]] = 2 * i;
		to[2 * i + 1] = g->eu[i]; nxt[2 * i + 1] = head[g->ev[i]]; head[g->ev[i]] = 2 * i + 1;
	}
	col = xmalloc((size_t)n * sizeof *col);
	stk = xmalloc((size_t)n * sizeof *stk);
	for (i = 0; i < n; i++)
		col[i] = -1;
	col[0] = 0;
	stk[top++] = 0;
	while (top && ok) {
		int v = stk[--top], a;
		for (a = head[v]; a >= 0; a = nxt[a]) {
			int w = to[a];
			if (col[w] < 0) {
				col[w] = 1 - col[v];
				stk[top++] = w;
			} else if (col[w] == col[v]) {
				ok = 0;			/* an odd cycle */
				break;
			}
		}
	}
	if (ok) {
		for (i = 0; i < n; i++)
			black += (col[i] == 1);
		*small = black < n - black ? black : n - black;
		if (out)
			for (i = 0; i < n; i++)
				out[i] = (unsigned char)(col[i] == 1);
	}
	free(head); free(nxt); free(to); free(col); free(stk);
	return ok;
}

/* Is the connected graph G a path?  Trees of maximum degree 2 and nothing else. */
static int is_path(const Graph *g)
{
	int n = g->n, i, ok = 1;
	int *deg;
	if (g->m != n - 1)
		return 0;
	deg = xcalloc((size_t)n, sizeof *deg);
	for (i = 0; i < g->m; i++) {
		deg[g->eu[i]]++;
		deg[g->ev[i]]++;
	}
	for (i = 0; i < n; i++)
		if (deg[i] > 2)
			ok = 0;
	free(deg);
	return ok;
}

/*
 * The hunt for a cut certificate, by the second lemma at the head of the file.
 * A candidate S is scored by (c - 1)/|S|, compared as a fraction so that the arithmetic stays integral, and ties go to the smaller S.
 * The candidates are, in this order and with no randomness anywhere, so that the route a graph takes is reproducible:
 *     both colour classes, when G is bipartite;
 *     the complement of a greedy maximal independent set, taking the vertices by ascending and then by descending degree;
 *     N(v), for every vertex v, which leaves v a component of its own.
 * The first family is exactly the bipartite bound, and the others reach graphs with odd cycles that it cannot see.
 */
typedef struct {
	const Graph *g;
	const int *head, *nxt, *to;
	unsigned char *seen;
	int *stack;
	unsigned char *best;
	int bsize, bcomp;
} CutHunt;

/* Count the components of G - S, with S given as vertex flags. */
static int cut_components(const Graph *g, const int *head, const int *nxt, const int *to,
			  const unsigned char *inS, unsigned char *seen, int *stack)
{
	int n = g->n, v, c = 0;
	memset(seen, 0, (size_t)n);
	for (v = 0; v < n; v++) {
		int top = 0;
		if (inS[v] || seen[v])
			continue;
		c++;
		seen[v] = 1;
		stack[top++] = v;
		while (top) {
			int u = stack[--top], a;
			for (a = head[u]; a >= 0; a = nxt[a]) {
				int w = to[a];
				if (!inS[w] && !seen[w]) {
					seen[w] = 1;
					stack[top++] = w;
				}
			}
		}
	}
	return c;
}

static void cut_offer(CutHunt *H, const unsigned char *cand)
{
	int n = H->g->n, i, sz = 0, c;
	long lhs, rhs;

	for (i = 0; i < n; i++)
		sz += cand[i];
	if (sz == 0 || sz == n)
		return;
	c = cut_components(H->g, H->head, H->nxt, H->to, cand, H->seen, H->stack);
	lhs = (long)(c - 1) * H->bsize;
	rhs = (long)(H->bcomp - 1) * sz;
	if (H->bsize == 0 || lhs > rhs || (lhs == rhs && sz < H->bsize)) {
		memcpy(H->best, cand, (size_t)n);
		H->bsize = sz;
		H->bcomp = c;
	}
}

/* Write the best candidate to best[], its size to *bsize and its component count to *bcomp.  *bsize is 0 when there is no candidate at all. */
static void cut_search(const Graph *g, int bip, const unsigned char *col,
		       unsigned char *best, int *bsize, int *bcomp)
{
	int n = g->n, m = g->m, i, j, v, a, pass;
	int *head, *nxt, *to, *deg, *ord, *stack;
	unsigned char *cand, *seen, *banned;
	CutHunt H;

	*bsize = 0;
	*bcomp = 0;
	memset(best, 0, (size_t)n);
	if (n < 3)
		return;

	head = xmalloc((size_t)n * sizeof *head);
	nxt = xmalloc((size_t)m * 2 * sizeof *nxt);
	to = xmalloc((size_t)m * 2 * sizeof *to);
	for (i = 0; i < n; i++)
		head[i] = -1;
	for (i = 0; i < m; i++) {
		to[2 * i] = g->ev[i]; nxt[2 * i] = head[g->eu[i]]; head[g->eu[i]] = 2 * i;
		to[2 * i + 1] = g->eu[i]; nxt[2 * i + 1] = head[g->ev[i]]; head[g->ev[i]] = 2 * i + 1;
	}
	deg = xcalloc((size_t)n, sizeof *deg);
	for (i = 0; i < m; i++) {
		deg[g->eu[i]]++;
		deg[g->ev[i]]++;
	}
	ord = xmalloc((size_t)n * sizeof *ord);
	for (i = 0; i < n; i++) {		/* by ascending degree, insertion sort: n is small */
		int key = i;
		for (j = i; j > 0 && deg[ord[j - 1]] > deg[key]; j--)
			ord[j] = ord[j - 1];
		ord[j] = key;
	}
	cand = xmalloc((size_t)n);
	seen = xmalloc((size_t)n);
	banned = xmalloc((size_t)n);
	stack = xmalloc((size_t)n * sizeof *stack);

	H.g = g; H.head = head; H.nxt = nxt; H.to = to;
	H.seen = seen; H.stack = stack; H.best = best;
	H.bsize = 0; H.bcomp = 0;

	if (bip) {
		for (i = 0; i < n; i++) cand[i] = !col[i];
		cut_offer(&H, cand);
		for (i = 0; i < n; i++) cand[i] = col[i];
		cut_offer(&H, cand);
	}
	for (pass = 0; pass < 2; pass++) {
		memset(banned, 0, (size_t)n);
		memset(cand, 1, (size_t)n);
		for (i = 0; i < n; i++) {
			v = pass ? ord[n - 1 - i] : ord[i];
			if (banned[v])
				continue;
			cand[v] = 0;		/* v joins the independent set, so it leaves S */
			banned[v] = 1;
			for (a = head[v]; a >= 0; a = nxt[a])
				banned[to[a]] = 1;
		}
		cut_offer(&H, cand);
	}
	for (v = 0; v < n; v++) {
		memset(cand, 0, (size_t)n);
		for (a = head[v]; a >= 0; a = nxt[a])
			cand[to[a]] = 1;
		cut_offer(&H, cand);
	}

	*bsize = H.bsize;
	*bcomp = H.bcomp;
	free(head); free(nxt); free(to); free(deg); free(ord);
	free(cand); free(seen); free(banned); free(stack);
}

/* ------------------------------------------------------------- the row store */

/*
 * The program is held here in integer form, independently of the solver, so that the exact pass can read the same rows the simplex was given.
 * Columns 0..m-1 are the x_e, and column m is the slack t.
 */
enum { ROW_LE, ROW_EQ };

typedef struct {
	int *idx, *cf;
	int len, rhs, sense;
} Row;

typedef struct {
	Row *r;
	int n, cap, ncol;
} Rows;

static void rows_free(Rows *R)
{
	int i;
	for (i = 0; i < R->n; i++) {
		free(R->r[i].idx);
		free(R->r[i].cf);
	}
	free(R->r);
}

static int rows_add(Rows *R, const int *idx, const int *cf, int len, int rhs, int sense)
{
	Row *row;
	if (R->n == R->cap) {
		R->cap = R->cap ? 2 * R->cap : 64;
		R->r = xrealloc(R->r, (size_t)R->cap * sizeof *R->r);
	}
	row = &R->r[R->n];
	row->idx = xmalloc((size_t)len * sizeof *row->idx);
	row->cf = xmalloc((size_t)len * sizeof *row->cf);
	memcpy(row->idx, idx, (size_t)len * sizeof *idx);
	memcpy(row->cf, cf, (size_t)len * sizeof *cf);
	row->len = len;
	row->rhs = rhs;
	row->sense = sense;
	return R->n++;
}

static double row_value(const Row *row, const double *z)
{
	double s = 0.0;
	int i;
	for (i = 0; i < row->len; i++)
		s += row->cf[i] * z[row->idx[i]];
	return s;
}

/*
 * Build the fixed part of the program.
 * Passing strict = 0 gives the RN degree rows d_v(x) <= 2, and strict = 1 gives the RP rows d_v(x) + t <= 2.
 */
static void rows_build(Rows *R, const Graph *g, int strict)
{
	int n = g->n, m = g->m, i, v, len;
	int *idx = xmalloc((size_t)(m + 2) * sizeof *idx);
	int *cf = xmalloc((size_t)(m + 2) * sizeof *cf);

	memset(R, 0, sizeof *R);
	R->ncol = m + 1;
	idx[0] = cf[0] = 0;		/* only the m = 0 case would read these, and that case never reaches here */

	len = 0;					/* x(E) = n - 1 */
	for (i = 0; i < m; i++) {
		idx[len] = i;
		cf[len++] = 1;
	}
	rows_add(R, idx, cf, len, n - 1, ROW_EQ);

	for (i = 0; i < m; i++) {			/* x_e >= t */
		idx[0] = i; cf[0] = -1;
		idx[1] = m; cf[1] = 1;
		rows_add(R, idx, cf, 2, 0, ROW_LE);
	}

	for (v = 0; v < n; v++) {			/* d_v(x) [+ t] <= 2 */
		len = 0;
		for (i = 0; i < m; i++)
			if (g->eu[i] == v || g->ev[i] == v) {
				idx[len] = i;
				cf[len++] = 1;
			}
		if (strict) { idx[len] = m; cf[len++] = 1; }
		rows_add(R, idx, cf, len, 2, ROW_LE);
	}

	for (i = 0; i < m; i++) {			/* x_e >= 0 */
		idx[0] = i; cf[0] = -1;
		rows_add(R, idx, cf, 1, 0, ROW_LE);
	}
	idx[0] = m; cf[0] = 1;				/* t <= 1 */
	rows_add(R, idx, cf, 1, 1, ROW_LE);
	idx[0] = m; cf[0] = -1;				/* t >= -n, which keeps it bounded */
	rows_add(R, idx, cf, 1, n, ROW_LE);

	free(idx);
	free(cf);
}

/* Add row `k` of the store to the solver. */
static void glp_push_row(glp_prob *lp, const Rows *R, int k)
{
	const Row *row = &R->r[k];
	int i, r = glp_add_rows(lp, 1);
	int *ind = xmalloc((size_t)(row->len + 1) * sizeof *ind);
	double *val = xmalloc((size_t)(row->len + 1) * sizeof *val);
	for (i = 0; i < row->len; i++) {
		ind[i + 1] = row->idx[i] + 1;
		val[i + 1] = row->cf[i];
	}
	if (row->sense == ROW_EQ)
		glp_set_row_bnds(lp, r, GLP_FX, row->rhs, row->rhs);
	else
		glp_set_row_bnds(lp, r, GLP_UP, 0.0, row->rhs);
	glp_set_mat_row(lp, r, row->len, ind, val);
	free(ind);
	free(val);
}

static glp_prob *glp_from_rows(const Rows *R, int n, int free_cols)
{
	glp_prob *lp = glp_create_prob();
	int m = R->ncol - 1, j;
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_cols(lp, R->ncol);
	/* Every bound is a row in its own right, so these column bounds only repeat what the program already says.
	   They are here because the simplex is quicker with them.
	   They have to go when a certificate is wanted: a variable resting on a column bound shows up as a reduced cost rather than as a row dual, and the exact pass picks the rows to rebuild the vertex from by row dual. */
	if (free_cols) {
		for (j = 0; j <= m; j++)
			glp_set_col_bnds(lp, j + 1, GLP_FR, 0.0, 0.0);
	} else {
		for (j = 0; j < m; j++)
			glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, n);
		glp_set_col_bnds(lp, m + 1, GLP_DB, -(double)n, 1.0);
	}
	glp_set_obj_coef(lp, m + 1, 1.0);
	for (j = 0; j < R->n; j++)
		glp_push_row(lp, R, j);
	return lp;
}

/* --------------------------------------------------------- Dinic maximum flow */

/*
 * Plain Dinic over doubles.
 * A second copy over GMP integers appears further down, for the exact rerun of the separation.
 * The two are deliberately kept apart rather than made generic, since between them they are 200 lines and only the arithmetic differs.
 */
typedef struct {
	int nv, na, cap;
	int *head, *next, *to, *level, *it;
	double *c;
} DNet;

static void dnet_init(DNet *N, int nv, int maxarcs)
{
	N->nv = nv;
	N->na = 0;
	N->cap = 2 * maxarcs;
	N->head = xmalloc((size_t)nv * sizeof *N->head);
	N->level = xmalloc((size_t)nv * sizeof *N->level);
	N->it = xmalloc((size_t)nv * sizeof *N->it);
	N->next = xmalloc((size_t)N->cap * sizeof *N->next);
	N->to = xmalloc((size_t)N->cap * sizeof *N->to);
	N->c = xmalloc((size_t)N->cap * sizeof *N->c);
}

static void dnet_free(DNet *N)
{
	free(N->head); free(N->next); free(N->to); free(N->c);
	free(N->level); free(N->it);
}

static void dnet_reset(DNet *N)
{
	int i;
	N->na = 0;
	for (i = 0; i < N->nv; i++)
		N->head[i] = -1;
}

static void dnet_arc(DNet *N, int u, int v, double cuv, double cvu)
{
	N->to[N->na] = v; N->c[N->na] = cuv; N->next[N->na] = N->head[u]; N->head[u] = N->na++;
	N->to[N->na] = u; N->c[N->na] = cvu; N->next[N->na] = N->head[v]; N->head[v] = N->na++;
}

static int dnet_bfs(DNet *N, int s, int t, int *q)
{
	int i, qh = 0, qt = 0;
	for (i = 0; i < N->nv; i++)
		N->level[i] = -1;
	N->level[s] = 0;
	q[qt++] = s;
	while (qh < qt) {
		int v = q[qh++], a;
		for (a = N->head[v]; a >= 0; a = N->next[a])
			if (N->c[a] > 1e-12 && N->level[N->to[a]] < 0) {
				N->level[N->to[a]] = N->level[v] + 1;
				q[qt++] = N->to[a];
			}
	}
	return N->level[t] >= 0;
}

static double dnet_dfs(DNet *N, int v, int t, double f)
{
	if (v == t)
		return f;
	for (; N->it[v] >= 0; N->it[v] = N->next[N->it[v]]) {
		int a = N->it[v], w = N->to[a];
		double d;
		if (N->c[a] <= 1e-12 || N->level[w] != N->level[v] + 1)
			continue;
		d = dnet_dfs(N, w, t, f < N->c[a] ? f : N->c[a]);
		if (d > 1e-12) {
			N->c[a] -= d;
			N->c[a ^ 1] += d;
			return d;
		}
	}
	N->level[v] = -1;
	return 0.0;
}

static double dnet_maxflow(DNet *N, int s, int t, int *q)
{
	double flow = 0.0;
	while (dnet_bfs(N, s, t, q)) {
		int i;
		double f;
		for (i = 0; i < N->nv; i++)
			N->it[i] = N->head[i];
		while ((f = dnet_dfs(N, s, t, 1e30)) > 1e-12)
			flow += f;
	}
	return flow;
}

/* Mark the source side of the minimum cut, as local vertex flags. */
static void dnet_cut(DNet *N, int s, int nloc, unsigned char *side, int *q)
{
	int i, qh = 0, qt = 0;
	for (i = 0; i < N->nv; i++)
		N->level[i] = 0;
	memset(side, 0, (size_t)nloc);
	N->level[s] = 1;
	q[qt++] = s;
	while (qh < qt) {
		int v = q[qh++], a;
		if (v < nloc)
			side[v] = 1;
		for (a = N->head[v]; a >= 0; a = N->next[a])
			if (N->c[a] > 1e-9 && !N->level[N->to[a]]) {
				N->level[N->to[a]] = 1;
				q[qt++] = N->to[a];
			}
	}
}

/* ---------------------------------------------- separating the rank inequalities */

/*
 * Find a violated rank inequality, if there is one.
 * The inequality is
 *
 *     x(E[S]) + (|S|-1) t <= |S|-1,   that is,   Phi(S) >= C,   C = 1 - t,
 *
 * for Phi(S) = C|S| - x(E[S]), and
 *
 *     2 Phi(S) = sum_{v in S} (2C - d_v(x)) + x(dS)
 *
 * is a cut function of S, so minimising it is a minimum s-t cut.
 * Singletons give Phi = C exactly and so never register, but S = V would, so the minimisation runs over proper nonempty S only.
 * Fix one vertex v0 and, for every other u, force u out of S and then into it; that covers all such S in 2(n - 1) flows.
 *
 * Append each set found to *cuts, as n flags.
 */
static int separate(const Graph *g, const double *x, double t,
		    unsigned char **cuts, int *ncuts, int *cap)
{
	int n = g->n, m = g->m, i, k, added = 0, s = n, snk = n + 1;
	double C = 1.0 - t, neg = 0.0, inf = 1.0;
	double *w = xmalloc((size_t)n * sizeof *w);
	int *q = xmalloc((size_t)(n + 2) * sizeof *q);
	unsigned char *flags = xmalloc((size_t)n);
	unsigned char *found = NULL;
	double *viol = NULL;
	int nfound = 0;
	DNet N;

	if (n < 3) { free(w); free(q); free(flags); return 0; }
	found = xmalloc((size_t)(2 * n) * (size_t)n);
	viol = xmalloc((size_t)(2 * n) * sizeof *viol);

	for (i = 0; i < n; i++)
		w[i] = 2.0 * C;
	for (i = 0; i < m; i++) {
		w[g->eu[i]] -= x[i];
		w[g->ev[i]] -= x[i];
		inf += 2.0 * fabs(x[i]);
	}
	for (i = 0; i < n; i++) {
		if (w[i] < 0.0) neg += -w[i];
		inf += fabs(w[i]);
	}

	dnet_init(&N, n + 2, m + 2 * n + 4);
	for (k = 0; k < 2 * (n - 1); k++) {
		int u = 1 + k / 2, in_, out_;
		double f;

		if (k % 2 == 0) { in_ = 0; out_ = u; }	/* v0 in S, u out */
		else            { in_ = u; out_ = 0; }	/* u in S, v0 out */

		dnet_reset(&N);
		for (i = 0; i < m; i++)
			if (x[i] > 1e-12)
				dnet_arc(&N, g->eu[i], g->ev[i], x[i], x[i]);
		for (i = 0; i < n; i++) {
			if (w[i] < -1e-15)     dnet_arc(&N, s, i, -w[i], 0.0);
			else if (w[i] > 1e-15) dnet_arc(&N, i, snk, w[i], 0.0);
		}
		dnet_arc(&N, s, in_, inf, 0.0);
		dnet_arc(&N, out_, snk, inf, 0.0);

		f = dnet_maxflow(&N, s, snk, q);
		if (f - neg > 2.0 * C - CUT_TOL)
			continue;

		dnet_cut(&N, s, n, flags, q);
		memcpy(found + (size_t)nfound * n, flags, (size_t)n);
		viol[nfound++] = 2.0 * C - (f - neg);
	}

	/* Every set found goes in, worst first.
	   Feeding the solver one row a round and letting it re-optimise between them is far slower: on P_15 x P_15 that is 270 rounds and half a minute, against 18 rounds and under two seconds. */
	for (;;) {
		int best = -1;
		for (i = 0; i < nfound; i++)
			if (viol[i] > CUT_TOL && (best < 0 || viol[i] > viol[best]))
				best = i;
		if (best < 0)
			break;
		viol[best] = 0.0;
		memcpy(flags, found + (size_t)best * n, (size_t)n);
		for (i = 0; i < *ncuts + added; i++)
			if (!memcmp(*cuts + (size_t)i * n, flags, (size_t)n))
				break;
		if (i < *ncuts + added)
			continue;
		if (*ncuts + added == *cap) {
			*cap = *cap ? 2 * *cap : 128;
			*cuts = xrealloc(*cuts, (size_t)*cap * (size_t)n);
		}
		memcpy(*cuts + (size_t)(*ncuts + added) * n, flags, (size_t)n);
		added++;
	}
	*ncuts += added;
	dnet_free(&N);
	free(w); free(q); free(flags);
	free(found); free(viol);
	return added;
}

/* Build the row for the set S. */
static void add_cut_row(Rows *R, const Graph *g, const unsigned char *flags)
{
	int n = g->n, m = R->ncol - 1, i, len = 0, size = 0;
	int *idx = xmalloc((size_t)(m + 2) * sizeof *idx);
	int *cf = xmalloc((size_t)(m + 2) * sizeof *cf);

	for (i = 0; i < n; i++)
		size += flags[i];
	for (i = 0; i < m; i++)
		if (flags[g->eu[i]] && flags[g->ev[i]]) {
			idx[len] = i;
			cf[len++] = 1;
		}
	idx[len] = m;
	cf[len++] = size - 1;
	rows_add(R, idx, cf, len, size - 1, ROW_LE);
	free(idx);
	free(cf);
}

/* ------------------------------------------------------------- one program */

typedef struct {
	Rows R;
	double *z;		/* the floating point optimum, R.ncol values */
	double *dual;		/* the row duals, for ranking the tight rows */
	double opt;		/* z[m], the slack */
	int st;			/* 0 optimal, 1 infeasible, -1 solver failure */
	int converged;		/* the separation ran out of violated rows */
	int ncuts, rounds;
} Program;

static void program_free(Program *P)
{
	rows_free(&P->R);
	free(P->z);
	free(P->dual);
}

/*
 * Solve by row generation.
 * Both programs are feasible, since t = -n does it, and bounded, since t <= 1, so an optimum always exists.
 * A failure here therefore belongs to the solver, and it is reported as such rather than as an answer.
 */
static void program_solve(Program *P, const Graph *g, int strict, int exact)
{
	glp_prob *lp;
	glp_smcp parm;
	unsigned char *cuts = NULL;
	int ncuts = 0, cap = 0, round, i;

	memset(P, 0, sizeof *P);
	P->st = -1;
	rows_build(&P->R, g, strict);
	P->z = xcalloc((size_t)P->R.ncol, sizeof *P->z);
	lp = glp_from_rows(&P->R, g->n, exact);

	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF;
	parm.presolve = GLP_OFF;

	for (round = 1; round <= MAX_ROUNDS; round++) {
		int added, before = ncuts, st;
		if (glp_simplex(lp, &parm) != 0)
			goto done;
		st = glp_get_status(lp);
		if (st == GLP_NOFEAS || st == GLP_INFEAS) {
			/* Every generated row is valid, so an infeasible relaxation means the true program is infeasible.
			   No point of P(G) then has all degrees at most 2, and G is not RN. */
			P->st = 1;
			goto done;
		}
		if (st != GLP_OPT)
			goto done;
		P->st = 0;
		for (i = 0; i < P->R.ncol; i++)
			P->z[i] = glp_get_col_prim(lp, i + 1);
		P->opt = P->z[P->R.ncol - 1];

		added = separate(g, P->z, P->opt, &cuts, &ncuts, &cap);
		if (!added) {
			P->converged = 1;
			break;
		}
		for (i = before; i < ncuts; i++) {
			add_cut_row(&P->R, g, cuts + (size_t)i * g->n);
			glp_push_row(lp, &P->R, P->R.n - 1);
		}
	}
done:
	if (P->st == 0) {
		P->dual = xcalloc((size_t)P->R.n, sizeof *P->dual);
		for (i = 0; i < P->R.n; i++)
			P->dual[i] = glp_get_row_dual(lp, i + 1);
	}
	P->ncuts = ncuts;
	P->rounds = round > MAX_ROUNDS ? MAX_ROUNDS : round;
	glp_delete_prob(lp);
	free(cuts);
}

/* ------------------------------------------------------- exact linear algebra */

/*
 * Solve the square integer system A z = b exactly.
 * Forward elimination is fraction free (Bareiss), so it stays in the integers and the entries stay bounded by minors of A.
 * The back substitution is over the rationals and costs only O(N^2).
 * Return 0 if A is singular.
 */
static int exact_solve(int N, const int *A, const int *b, mpq_t *z)
{
	mpz_t *M = xmalloc((size_t)N * (size_t)(N + 1) * sizeof *M);
	mpz_t prev, tmp;
	mpq_t acc, qt;
	int i, j, k, ok = 1;
#define MM(r, c) M[(size_t)(r) * (N + 1) + (c)]

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			mpz_init_set_si(MM(i, j), A[(size_t)i * N + j]);
		mpz_init_set_si(MM(i, N), b[i]);
	}
	mpz_init_set_ui(prev, 1);
	mpz_init(tmp);

	for (k = 0; k < N && ok; k++) {
		int p = -1;
		for (i = k; i < N; i++)
			if (mpz_sgn(MM(i, k)) != 0) { p = i; break; }
		if (p < 0) { ok = 0; break; }
		if (p != k)
			for (j = k; j <= N; j++)
				mpz_swap(MM(k, j), MM(p, j));
		for (i = k + 1; i < N; i++) {
			if (mpz_sgn(MM(i, k)) == 0) {
				for (j = k + 1; j <= N; j++) {
					mpz_mul(tmp, MM(i, j), MM(k, k));
					mpz_divexact(MM(i, j), tmp, prev);
				}
				continue;
			}
			for (j = k + 1; j <= N; j++) {
				mpz_mul(tmp, MM(i, j), MM(k, k));
				mpz_submul(tmp, MM(i, k), MM(k, j));
				mpz_divexact(MM(i, j), tmp, prev);
			}
			mpz_set_ui(MM(i, k), 0);
		}
		mpz_set(prev, MM(k, k));
	}

	if (ok) {
		mpq_init(acc);
		mpq_init(qt);
		for (i = N - 1; i >= 0; i--) {
			mpq_set_z(acc, MM(i, N));
			for (j = i + 1; j < N; j++) {
				if (mpz_sgn(MM(i, j)) == 0 || mpq_sgn(z[j]) == 0)
					continue;
				mpq_set_z(qt, MM(i, j));
				mpq_mul(qt, qt, z[j]);
				mpq_sub(acc, acc, qt);
			}
			mpq_set_z(qt, MM(i, i));
			mpq_div(z[i], acc, qt);
		}
		mpq_clear(acc);
		mpq_clear(qt);
	}

	for (i = 0; i < N * (N + 1); i++)
		mpz_clear(M[i]);
	free(M);
	mpz_clear(prev);
	mpz_clear(tmp);
	return ok;
#undef MM
}

/* --------------------------------------------- Dinic over GMP integers */

typedef struct {
	int nv, na, cap;
	int *head, *next, *to, *level, *it;
	mpz_t *c;
} ZNet;

static void znet_init(ZNet *N, int nv, int maxarcs)
{
	int i;
	N->nv = nv;
	N->na = 0;
	N->cap = 2 * maxarcs;
	N->head = xmalloc((size_t)nv * sizeof *N->head);
	N->level = xmalloc((size_t)nv * sizeof *N->level);
	N->it = xmalloc((size_t)nv * sizeof *N->it);
	N->next = xmalloc((size_t)N->cap * sizeof *N->next);
	N->to = xmalloc((size_t)N->cap * sizeof *N->to);
	N->c = xmalloc((size_t)N->cap * sizeof *N->c);
	for (i = 0; i < N->cap; i++)
		mpz_init(N->c[i]);
}

static void znet_free(ZNet *N)
{
	int i;
	for (i = 0; i < N->cap; i++)
		mpz_clear(N->c[i]);
	free(N->c); free(N->head); free(N->next); free(N->to);
	free(N->level); free(N->it);
}

static void znet_reset(ZNet *N)
{
	int i;
	N->na = 0;
	for (i = 0; i < N->nv; i++)
		N->head[i] = -1;
}

static void znet_arc(ZNet *N, int u, int v, const mpz_t cuv, const mpz_t cvu)
{
	mpz_set(N->c[N->na], cuv);
	N->to[N->na] = v; N->next[N->na] = N->head[u]; N->head[u] = N->na++;
	mpz_set(N->c[N->na], cvu);
	N->to[N->na] = u; N->next[N->na] = N->head[v]; N->head[v] = N->na++;
}

static int znet_bfs(ZNet *N, int s, int t, int *q)
{
	int i, qh = 0, qt = 0;
	for (i = 0; i < N->nv; i++)
		N->level[i] = -1;
	N->level[s] = 0;
	q[qt++] = s;
	while (qh < qt) {
		int v = q[qh++], a;
		for (a = N->head[v]; a >= 0; a = N->next[a])
			if (mpz_sgn(N->c[a]) > 0 && N->level[N->to[a]] < 0) {
				N->level[N->to[a]] = N->level[v] + 1;
				q[qt++] = N->to[a];
			}
	}
	return N->level[t] >= 0;
}

static void znet_dfs(ZNet *N, int v, int t, const mpz_t f, mpz_t out)
{
	if (v == t) { mpz_set(out, f); return; }
	for (; N->it[v] >= 0; N->it[v] = N->next[N->it[v]]) {
		int a = N->it[v], w = N->to[a];
		mpz_t lim, d;
		if (mpz_sgn(N->c[a]) <= 0 || N->level[w] != N->level[v] + 1)
			continue;
		mpz_init(lim);
		mpz_init(d);
		mpz_set(lim, mpz_cmp(f, N->c[a]) < 0 ? f : N->c[a]);
		znet_dfs(N, w, t, lim, d);
		if (mpz_sgn(d) > 0) {
			mpz_sub(N->c[a], N->c[a], d);
			mpz_add(N->c[a ^ 1], N->c[a ^ 1], d);
			mpz_set(out, d);
			mpz_clear(lim);
			mpz_clear(d);
			return;
		}
		mpz_clear(lim);
		mpz_clear(d);
	}
	N->level[v] = -1;
	mpz_set_ui(out, 0);
}

static void znet_maxflow(ZNet *N, int s, int t, int *q, const mpz_t big, mpz_t flow)
{
	mpz_t f;
	mpz_init(f);
	mpz_set_ui(flow, 0);
	while (znet_bfs(N, s, t, q)) {
		int i;
		for (i = 0; i < N->nv; i++)
			N->it[i] = N->head[i];
		for (;;) {
			znet_dfs(N, s, t, big, f);
			if (mpz_sgn(f) <= 0)
				break;
			mpz_add(flow, flow, f);
		}
	}
	mpz_clear(f);
}

/* ------------------------------------------------------------- certification */

typedef struct {
	int ncol;
	int primal;		/* the exact vertex was rebuilt and is feasible */
	int dual;		/* the exact dual is sign feasible */
	mpq_t *z;		/* the exact vertex, valid when primal */
	mpq_t lb;		/* z[t]: a certified lower bound on the optimum */
	mpq_t ub;		/* b.y: a certified upper bound on it */
} Cert;

static void cert_init(Cert *C, int ncol)
{
	int i;
	memset(C, 0, sizeof *C);
	C->ncol = ncol;
	C->z = xmalloc((size_t)ncol * sizeof *C->z);
	for (i = 0; i < ncol; i++)
		mpq_init(C->z[i]);
	mpq_init(C->lb);
	mpq_init(C->ub);
}

static void cert_free(Cert *C)
{
	int i;
	for (i = 0; i < C->ncol; i++)
		mpq_clear(C->z[i]);
	free(C->z);
	mpq_clear(C->lb);
	mpq_clear(C->ub);
}

/*
 * Collect the rows tight at the floating point optimum, thinned down to a square nonsingular system.
 * Independence is tested in double precision, so the selection is only a guess at which rows to pick.
 * The exact solve below rejects a singular choice, and everything the choice leads to is checked exactly afterwards.
 */
typedef struct { double key; int idx; } Ranked;

static int ranked_cmp(const void *a, const void *b)
{
	const Ranked *x = a, *y = b;
	if (x->key != y->key)
		return x->key > y->key ? -1 : 1;
	return x->idx < y->idx ? -1 : (x->idx > y->idx);
}

static int select_tight(const Rows *R, const double *z, const double *dual, int *sel)
{
	int ncol = R->ncol, i, c, k = 0;
	double **piv = xcalloc((size_t)ncol, sizeof *piv);
	double *tmp = xmalloc((size_t)ncol * sizeof *tmp);
	Ranked *ord = xmalloc((size_t)R->n * sizeof *ord);

	for (i = 0; i < R->n; i++) {
		ord[i].idx = i;
		ord[i].key = dual ? fabs(dual[i]) : 0.0;
	}
	qsort(ord, (size_t)R->n, sizeof *ord, ranked_cmp);

	for (c = 0; c < R->n && k < ncol; c++) {
		const Row *row;
		i = ord[c].idx;
		row = &R->r[i];
		double val = row_value(row, z);
		double scale = abs(row->rhs) > 1 ? (double)abs(row->rhs) : 1.0;
		int j, best = -1;

		if (fabs(val - row->rhs) > FEAS_TOL * scale)
			continue;
		memset(tmp, 0, (size_t)ncol * sizeof *tmp);
		for (j = 0; j < row->len; j++)
			tmp[row->idx[j]] += row->cf[j];
		for (j = 0; j < ncol; j++) {
			if (fabs(tmp[j]) < 1e-9)
				continue;
			if (piv[j]) {
				double fq = tmp[j];
				int q;
				for (q = j; q < ncol; q++)
					tmp[q] -= fq * piv[j][q];
			} else {
				best = j;
				break;
			}
		}
		if (best < 0)
			continue;			/* dependent on what we have */
		{
			double f = tmp[best];
			for (j = best; j < ncol; j++)
				tmp[j] /= f;
		}
		piv[best] = xmalloc((size_t)ncol * sizeof *piv[best]);
		memcpy(piv[best], tmp, (size_t)ncol * sizeof *tmp);
		sel[k++] = i;
	}
	for (c = 0; c < ncol; c++)
		free(piv[c]);
	free(piv);
	free(tmp);
	free(ord);
	return k;
}

/* Does the exact point satisfy every row the program carries? */
static int exact_rows_ok(const Rows *R, mpq_t *z)
{
	mpq_t acc, term, rhs;
	int i, j, ok = 1;
	mpq_init(acc);
	mpq_init(term);
	mpq_init(rhs);
	for (i = 0; i < R->n && ok; i++) {
		const Row *row = &R->r[i];
		int cmp;
		mpq_set_ui(acc, 0, 1);
		for (j = 0; j < row->len; j++) {
			mpq_set_si(term, row->cf[j], 1);
			mpq_mul(term, term, z[row->idx[j]]);
			mpq_add(acc, acc, term);
		}
		mpq_set_si(rhs, row->rhs, 1);
		cmp = mpq_cmp(acc, rhs);
		ok = (row->sense == ROW_EQ) ? (cmp == 0) : (cmp <= 0);
	}
	mpq_clear(acc);
	mpq_clear(term);
	mpq_clear(rhs);
	return ok;
}

/*
 * Rerun exactly the rank inequalities the program does not carry.
 * Clearing denominators by D turns the cut function into integers: with X_e = D x_e, T = D t and C = D - T,
 *
 *     2 (C|S| - X(E_B[S])) = sum_{v in S} (2C - D_v(X)) + X(dS)  >=  2C
 *
 * has to hold for every proper nonempty S of vertices.
 */
static int exact_separate_ok(const Graph *g, mpq_t *z)
{
	int m = g->m, n = g->n, i, k, s = n, snk = n + 1, ok = 1;
	mpz_t D, C, twoC, big, neg, flow, t0;
	mpz_t *X = xmalloc((size_t)(m + 1) * sizeof *X);
	mpz_t *w = xmalloc((size_t)n * sizeof *w);
	int *q = xmalloc((size_t)(n + 2) * sizeof *q);
	ZNet N;

	mpz_init_set_ui(D, 1);
	for (i = 0; i <= m; i++) {
		mpz_init(X[i]);
		mpz_lcm(D, D, mpq_denref(z[i]));
	}
	for (i = 0; i <= m; i++) {		/* X_i = D z_i, exactly */
		mpz_mul(X[i], D, mpq_numref(z[i]));
		mpz_divexact(X[i], X[i], mpq_denref(z[i]));
	}
	mpz_init(C);
	mpz_sub(C, D, X[m]);
	mpz_init(twoC);
	mpz_mul_ui(twoC, C, 2);
	mpz_init(big);
	mpz_init(neg);
	mpz_init(flow);
	mpz_init(t0);

	for (i = 0; i < n; i++) {
		mpz_init(w[i]);
		mpz_set(w[i], twoC);
	}
	mpz_set_ui(big, 1);
	for (i = 0; i < m; i++) {
		mpz_sub(w[g->eu[i]], w[g->eu[i]], X[i]);
		mpz_sub(w[g->ev[i]], w[g->ev[i]], X[i]);
		mpz_abs(t0, X[i]);
		mpz_addmul_ui(big, t0, 2);
	}
	mpz_set_ui(neg, 0);
	for (i = 0; i < n; i++) {
		mpz_abs(t0, w[i]);
		mpz_add(big, big, t0);
		if (mpz_sgn(w[i]) < 0)
			mpz_add(neg, neg, t0);
	}
	mpz_add(big, big, twoC);

	if (n >= 3) {
		znet_init(&N, n + 2, m + 2 * n + 4);
		for (k = 0; k < 2 * (n - 1) && ok; k++) {
			int u = 1 + k / 2;
			int in_ = (k % 2 == 0) ? 0 : u;
			int out_ = (k % 2 == 0) ? u : 0;

			znet_reset(&N);
			for (i = 0; i < m; i++)
				if (mpz_sgn(X[i]) > 0)
					znet_arc(&N, g->eu[i], g->ev[i], X[i], X[i]);
			mpz_set_ui(t0, 0);
			for (i = 0; i < n; i++) {
				if (mpz_sgn(w[i]) < 0) {
					mpz_neg(t0, w[i]);
					znet_arc(&N, s, i, t0, t0);
					mpz_set_ui(N.c[N.na - 1], 0);
				} else if (mpz_sgn(w[i]) > 0) {
					znet_arc(&N, i, snk, w[i], w[i]);
					mpz_set_ui(N.c[N.na - 1], 0);
				}
			}
			znet_arc(&N, s, in_, big, big);
			mpz_set_ui(N.c[N.na - 1], 0);
			znet_arc(&N, out_, snk, big, big);
			mpz_set_ui(N.c[N.na - 1], 0);

			znet_maxflow(&N, s, snk, q, big, flow);
			mpz_sub(flow, flow, neg);	/* = 2 min Phi over these S */
			if (mpz_cmp(flow, twoC) < 0)
				ok = 0;
		}
		znet_free(&N);
	}

	for (i = 0; i < n; i++)
		mpz_clear(w[i]);
	free(w);
	free(q);
	for (i = 0; i <= m; i++)
		mpz_clear(X[i]);
	free(X);
	mpz_clear(D); mpz_clear(C); mpz_clear(twoC); mpz_clear(big);
	mpz_clear(neg); mpz_clear(flow); mpz_clear(t0);
	return ok;
}

/*
 * Rebuild the vertex and its dual exactly from the same square system, and check both.
 * A verified vertex with t > 0 proves the property, and a verified dual with b.y <= 0 disproves it.
 * Since y is supported on the selected rows, A_sel^T y = c already says A^T y = c.
 * For any feasible x it follows that t = c.x = y.(A_sel x) <= y.b_sel, and that is the bound.
 */
static void certify(Cert *C, const Graph *g, const Program *P)
{
	int ncol = P->R.ncol, i, j;
	int *sel = xmalloc((size_t)ncol * sizeof *sel);
	int *A, *bb, nsel;
	mpq_t *y, term;

	cert_init(C, ncol);
	nsel = select_tight(&P->R, P->z, P->dual, sel);
	if (nsel < ncol) { free(sel); return; }

	A = xcalloc((size_t)ncol * (size_t)ncol, sizeof *A);
	bb = xmalloc((size_t)ncol * sizeof *bb);
	for (i = 0; i < ncol; i++) {
		const Row *row = &P->R.r[sel[i]];
		for (j = 0; j < row->len; j++)
			A[(size_t)i * ncol + row->idx[j]] += row->cf[j];
		bb[i] = row->rhs;
	}

	if (exact_solve(ncol, A, bb, C->z)) {
		if (exact_rows_ok(&P->R, C->z) && exact_separate_ok(g, C->z)) {
			C->primal = 1;
			mpq_set(C->lb, C->z[ncol - 1]);
		}
	}

	/* the transposed system, for the dual */
	y = xmalloc((size_t)ncol * sizeof *y);
	for (i = 0; i < ncol; i++)
		mpq_init(y[i]);
	{
		int *At = xmalloc((size_t)ncol * (size_t)ncol * sizeof *At);
		int *c = xcalloc((size_t)ncol, sizeof *c);
		for (i = 0; i < ncol; i++)
			for (j = 0; j < ncol; j++)
				At[(size_t)i * ncol + j] = A[(size_t)j * ncol + i];
		c[ncol - 1] = 1;
		if (exact_solve(ncol, At, c, y)) {
			int ok = 1;
			for (i = 0; i < ncol; i++)
				if (P->R.r[sel[i]].sense == ROW_LE && mpq_sgn(y[i]) < 0)
					ok = 0;
			if (ok) {
				mpq_init(term);
				mpq_set_ui(C->ub, 0, 1);
				for (i = 0; i < ncol; i++) {
					mpq_set_si(term, P->R.r[sel[i]].rhs, 1);
					mpq_mul(term, term, y[i]);
					mpq_add(C->ub, C->ub, term);
				}
				mpq_clear(term);
				C->dual = 1;
			}
		}
		free(At);
		free(c);
	}
	for (i = 0; i < ncol; i++)
		mpq_clear(y[i]);
	free(y);
	free(A);
	free(bb);
	free(sel);
}

/* -------------------------------------------------------------------- driver */

typedef struct {
	int witness, echo, g6, jsonl, exact;
} Opts;

/* Did rational arithmetic settle the sign of this program's optimum? */
static int cert_settled(const Program *P, const Cert *C)
{
	if (P->st != 0)
		return 0;
	if (C->primal && mpq_sgn(C->lb) > 0)
		return 1;
	return C->dual && mpq_sgn(C->ub) <= 0;
}

/* Return +1 if the optimum is positive, 0 if it is not, and -1 if the sign is unsettled. */
static int verdict(const Program *P, const Cert *C, int exact)
{
	if (P->st == 1)
		return 0;			/* infeasible: no such point at all */
	if (P->st != 0)
		return -1;
	if (exact) {
		if (C->primal && mpq_sgn(C->lb) > 0)
			return 1;
		if (C->dual && mpq_sgn(C->ub) <= 0)
			return 0;
		return -1;
	}
	return P->opt > OPT_TOL;
}

static void emit_vertex(FILE *f, const Graph *g, int v)
{
	if (g->vname)
		jsquote(f, g->vname[v]);
	else
		fprintf(f, "%d", v);
}

static void emit_mpq(FILE *f, const mpq_t q)
{
	char *s = mpq_get_str(NULL, 10, q);
	jsquote(f, s);
	free(s);
}

/* Report the witness in doubles, with its degrees and curvatures. */
static void emit_witness(FILE *f, const Graph *g, const double *x, int converged)
{
	int n = g->n, m = g->m, i, v, first = 1;
	double *d = xcalloc((size_t)n, sizeof *d);

	for (i = 0; i < m; i++) {
		d[g->eu[i]] += x[i];
		d[g->ev[i]] += x[i];
	}
	fputs("{\"point\":[", f);
	for (i = 0; i < m; i++) { if (i) fputs(",", f); jsnum(f, x[i]); }
	fputs("],\"degree\":[", f);
	for (v = 0; v < n; v++) { if (v) fputs(",", f); jsnum(f, d[v]); }
	fputs("],\"curvature\":[", f);
	for (v = 0; v < n; v++) { if (v) fputs(",", f); jsnum(f, 1.0 - d[v] / 2.0); }
	fputs("],\"zero_curvature\":[", f);
	for (v = 0; v < n; v++)
		if (fabs(1.0 - d[v] / 2.0) <= FEAS_TOL) {
			if (!first) fputs(",", f);
			first = 0;
			emit_vertex(f, g, v);
		}
	fprintf(f, "],\"verified\":%s}", converged ? "true" : "false");
	free(d);
}

/* Report the same in exact rationals. */
static void emit_witness_exact(FILE *f, const Graph *g, mpq_t *z)
{
	int n = g->n, m = g->m, i, v, first = 1;
	mpq_t *d = xmalloc((size_t)n * sizeof *d), half, one, p;

	for (v = 0; v < n; v++)
		mpq_init(d[v]);
	for (i = 0; i < m; i++) {
		mpq_add(d[g->eu[i]], d[g->eu[i]], z[i]);
		mpq_add(d[g->ev[i]], d[g->ev[i]], z[i]);
	}
	mpq_init(half); mpq_set_ui(half, 1, 2);
	mpq_init(one);  mpq_set_ui(one, 1, 1);
	mpq_init(p);

	fputs("{\"point\":[", f);
	for (i = 0; i < m; i++) { if (i) fputs(",", f); emit_mpq(f, z[i]); }
	fputs("],\"degree\":[", f);
	for (v = 0; v < n; v++) { if (v) fputs(",", f); emit_mpq(f, d[v]); }
	fputs("],\"curvature\":[", f);
	for (v = 0; v < n; v++) {
		if (v) fputs(",", f);
		mpq_mul(p, d[v], half);
		mpq_sub(p, one, p);
		emit_mpq(f, p);
	}
	fputs("],\"zero_curvature\":[", f);
	for (v = 0; v < n; v++) {
		mpq_mul(p, d[v], half);
		mpq_sub(p, one, p);
		if (mpq_sgn(p) == 0) {
			if (!first) fputs(",", f);
			first = 0;
			emit_vertex(f, g, v);
		}
	}
	fputs("]}", f);

	for (v = 0; v < n; v++)
		mpq_clear(d[v]);
	free(d);
	mpq_clear(half); mpq_clear(one); mpq_clear(p);
}

/*
 * Answer a graph that is not 2-connected, by the lemma at the head of this file: it is RN exactly when it is a path.
 * On a path every edge is a bridge, so P(G) = {1} and the witness is x = 1: the degrees are 2 at the internal vertices and 1 at the ends.
 * The RP program would then report t = min(1, 2 - max_v d_v), which is 1 for n <= 2 and 0 beyond, and the RN program t = 1.
 * Nothing else here is RN, and the certificate for that is the cut vertex itself.
 */
static void emit_shortcut(FILE *f, const Graph *g, int cut, const Opts *o)
{
	int n = g->n, m = g->m, i, path = is_path(g), rp = path && n <= 2;
	const char *cls = !path ? "not RN" : (rp ? "RP" : "SRN");

	fputs(",\"two_connected\":false,\"shortcut\":", f);
	jsquote(f, path ? "path" : "cut-vertex");
	if (!path && cut >= 0) {
		fputs(",\"cut_vertex\":", f);
		emit_vertex(f, g, cut);
	}
	fputs(",\"class\":", f);
	jsquote(f, cls);
	fprintf(f, ",\"rn\":%s,\"rp\":%s,\"srn\":%s",
		path ? "true" : "false", rp ? "true" : "false",
		(path && !rp) ? "true" : "false");
	fprintf(f, ",\"delta\":%d,\"lambda\":%d", rp ? 1 : 0, (path && !rp) ? 1 : 0);
	fputs(",\"rank_cuts\":0,\"rounds\":0", f);

	if (o->exact) {
		fputs(",\"exact\":{\"status\":\"certified\"", f);
		if (path) {
			/* Whichever of the two programs settles the class, its optimum here is 1. */
			fputs(",\"bound_lower\":\"1\",\"bound_upper\":\"1\",\"optimum\":\"1\"", f);
			if (o->witness) {
				mpq_t *z = xmalloc((size_t)(m + 1) * sizeof *z);
				for (i = 0; i <= m; i++) {
					mpq_init(z[i]);
					mpq_set_ui(z[i], 1, 1);
				}
				fputs(",\"witness\":", f);
				emit_witness_exact(f, g, z);
				for (i = 0; i <= m; i++)
					mpq_clear(z[i]);
				free(z);
			}
		} else {
			fputs(",\"reason\":", f);
			jsquote(f, "a cut vertex, and G is not a path: see the lemma at the head of rn.c");
		}
		fputs("}", f);
	}

	if (o->witness) {
		fputs(",\"witness\":", f);
		if (path) {
			double *x = xmalloc((size_t)(m + 1) * sizeof *x);
			for (i = 0; i <= m; i++)
				x[i] = 1.0;
			emit_witness(f, g, x, 1);
			free(x);
		} else {
			fputs("null", f);
		}
	}
}

/* The cut certificate itself, which a reader can check by hand. */
static void emit_cut(FILE *f, const Graph *g, const unsigned char *S, int size, int comps)
{
	int v, first = 1;
	fputs(",\"cut_set\":[", f);
	for (v = 0; v < g->n; v++)
		if (S[v]) {
			if (!first) fputs(",", f);
			first = 0;
			emit_vertex(f, g, v);
		}
	fprintf(f, "],\"cut_size\":%d,\"cut_components\":%d", size, comps);
}

/*
 * Answer a graph carrying a cut with c >= |S| + 2, by the second lemma at the head of this file.
 * There max_{v in S} d_v(x) >= 1 + (c-1)/|S| > 2 on all of P(G), so no point of it has every degree at most 2.
 */
static void emit_cut_gap(FILE *f, const Graph *g, const unsigned char *S, int size,
			 int comps, const Opts *o)
{
	fputs(",\"two_connected\":true,\"shortcut\":\"cut-gap\"", f);
	emit_cut(f, g, S, size, comps);
	fputs(",\"class\":\"not RN\",\"rn\":false,\"rp\":false,\"srn\":false"
	      ",\"delta\":null,\"lambda\":null,\"rank_cuts\":0,\"rounds\":0", f);
	if (o->exact) {
		char buf[192];
		snprintf(buf, sizeof buf,
			 "a cut of %d vertices leaves %d components, so max_v d_v(x) >= 1 + %d/%d > 2 on all of P(G)",
			 size, comps, comps - 1, size);
		fputs(",\"exact\":{\"status\":\"certified\",\"reason\":", f);
		jsquote(f, buf);
		fputs("}", f);
	}
	if (o->witness)
		fputs(",\"witness\":null", f);
}

static void run_one(FILE *f, Graph *g, const Opts *o)
{
	int n = g->n, m = g->m, i, v, cut, bip, small, skip_rp = 0;
	int csize = 0, ccomp = 0;
	unsigned char *col = NULL, *cutS = NULL;
	Program rpP, rnP;
	Cert rpC, rnC;
	int rp_v, rn_v = 0, rp = 0, rn = 0, ranrn = 0, cert_ok = 1, theta_empty = 0;
	const Program *W = NULL;		/* the program holding the witness */
	const Cert *WC = NULL;
	const char *cls;
	double t0, t1;
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);
	t0 = ts.tv_sec + 1e-9 * ts.tv_nsec;

	fputs("{", f);
	if (g->name) { fputs("\"name\":", f); jsquote(f, g->name); fputs(",", f); }
	fprintf(f, "\"n\":%d,\"m\":%d", n, m);

	if (o->echo) {
		fputs(",\"vertices\":[", f);
		for (v = 0; v < n; v++) { if (v) fputs(",", f); emit_vertex(f, g, v); }
		fputs("],\"edges\":[", f);
		for (i = 0; i < m; i++) {
			if (i) fputs(",", f);
			fputs("[", f);
			emit_vertex(f, g, g->eu[i]);
			fputs(",", f);
			emit_vertex(f, g, g->ev[i]);
			fputs("]", f);
		}
		fputs("]", f);
	}
	if (o->g6) {
		char *s = g6_encode(g);
		if (s) { fputs(",\"graph6\":", f); jsquote(f, s); free(s); }
	}
	if (g->note) { fputs(",\"note\":", f); jsquote(f, g->note); }

	if (n == 0) {
		fputs(",\"connected\":false,\"class\":null,\"error\":\"empty graph\"}", f);
		return;
	}
	if (!connected(g)) {
		fputs(",\"connected\":false,\"class\":null,"
		      "\"error\":\"resistance curvature is defined on connected graphs\"}", f);
		return;
	}
	fputs(",\"connected\":true", f);

	col = xmalloc((size_t)n);
	cutS = xmalloc((size_t)n);
	bip = bipartition(g, &small, col);
	fprintf(f, ",\"bipartite\":%s", bip ? "true" : "false");
	if (bip)
		fprintf(f, ",\"parts\":[%d,%d]", small, n - small);

	/* By the first lemma at the head of this file a graph that is not 2-connected needs no program at all. */
	if (!two_connected(g, &cut)) {
		emit_shortcut(f, g, cut, o);
		clock_gettime(CLOCK_MONOTONIC, &ts);
		t1 = ts.tv_sec + 1e-9 * ts.tv_nsec;
		fputs(",\"seconds\":", f);
		fprintf(f, "%.6f", t1 - t0);
		fputs("}", f);
		free(col); free(cutS);
		return;
	}
	/* And by the second, a cut leaving |S| + 2 components or more says not RN outright, while one leaving exactly |S| + 1 rules out RP and leaves only the RN program to run. */
	cut_search(g, bip, col, cutS, &csize, &ccomp);
	if (ccomp >= csize + 2) {
		emit_cut_gap(f, g, cutS, csize, ccomp, o);
		clock_gettime(CLOCK_MONOTONIC, &ts);
		t1 = ts.tv_sec + 1e-9 * ts.tv_nsec;
		fputs(",\"seconds\":", f);
		fprintf(f, "%.6f", t1 - t0);
		fputs("}", f);
		free(col); free(cutS);
		return;
	}
	skip_rp = (ccomp == csize + 1);
	fputs(",\"two_connected\":true,\"shortcut\":", f);
	if (skip_rp)
		jsquote(f, "cut-tight");
	else
		fputs("null", f);
	if (skip_rp)
		emit_cut(f, g, cutS, csize, ccomp);

	if (skip_rp) {
		/* max_{v in S} d_v(x) >= 2 throughout P(G), so the RP program would report an optimum of at most 0 and is not run. */
		memset(&rpP, 0, sizeof rpP);
		rpP.st = -1;
		cert_init(&rpC, m + 1);
		rp = rn = rp_v = 0;
	} else {
	program_solve(&rpP, g, 1, o->exact);
	if (o->exact && rpP.st == 0)
		certify(&rpC, g, &rpP);
	else
		cert_init(&rpC, rpP.R.ncol);
	rp_v = verdict(&rpP, &rpC, o->exact);
	if (rp_v < 0)
		rp_v = (rpP.st == 0 && rpP.opt > OPT_TOL);
	if (o->exact && !cert_settled(&rpP, &rpC))
		cert_ok = 0;
	rp = rn = rp_v;
	if (rp) { W = &rpP; WC = &rpC; }
	if (!rp) {
		if (o->exact && rpC.dual && mpq_sgn(rpC.ub) < 0)
			theta_empty = 1;
		else if (rpP.st == 0 && rpP.opt < -OPT_TOL)
			theta_empty = 1;
	}
	}

	if (!rp && !theta_empty) {
		ranrn = 1;
		program_solve(&rnP, g, 0, o->exact);
		if (o->exact && rnP.st == 0)
			certify(&rnC, g, &rnP);
		else
			cert_init(&rnC, rnP.R.ncol);
		rn_v = verdict(&rnP, &rnC, o->exact);
		if (rn_v < 0)
			rn_v = (rnP.st == 0 && rnP.opt > OPT_TOL);
		if (o->exact && !cert_settled(&rnP, &rnC))
			cert_ok = 0;
		rn = rn_v;
		if (rn) { W = &rnP; WC = &rnC; }
	}

	cls = rp ? "RP" : (rn ? "SRN" : "not RN");
	fputs(",\"class\":", f);
	jsquote(f, cls);
	fprintf(f, ",\"rn\":%s,\"rp\":%s,\"srn\":%s",
		rn ? "true" : "false", rp ? "true" : "false",
		(rn && !rp) ? "true" : "false");
	fputs(",\"delta\":", f);			/* null when the RP program was not run at all */
	if (skip_rp)
		fputs("null", f);
	else
		jsnum(f, rpP.st == 0 ? rpP.opt : 0.0);
	fputs(",\"lambda\":", f);
	jsnum(f, (ranrn && rnP.st == 0) ? rnP.opt : 0.0);
	fprintf(f, ",\"rank_cuts\":%d,\"rounds\":%d",
		rpP.ncuts + (ranrn ? rnP.ncuts : 0),
		rpP.rounds + (ranrn ? rnP.rounds : 0));

	if (o->exact) {
		const Cert *C = ranrn ? &rnC : &rpC;
		fputs(",\"exact\":{\"status\":", f);
		jsquote(f, cert_ok ? "certified" : "unresolved");
		if (C->primal) { fputs(",\"bound_lower\":", f); emit_mpq(f, C->lb); }
		if (C->dual)   { fputs(",\"bound_upper\":", f); emit_mpq(f, C->ub); }
		if (C->primal && C->dual && mpq_equal(C->lb, C->ub)) {
			fputs(",\"optimum\":", f);
			emit_mpq(f, C->lb);
		}
		if (o->witness && W && WC->primal) {
			fputs(",\"witness\":", f);
			emit_witness_exact(f, g, WC->z);
		}
		fputs("}", f);
	}

	if (o->witness) {
		fputs(",\"witness\":", f);
		if (W)
			emit_witness(f, g, W->z, W->converged);
		else
			fputs("null", f);
	}

	clock_gettime(CLOCK_MONOTONIC, &ts);
	t1 = ts.tv_sec + 1e-9 * ts.tv_nsec;
	fputs(",\"seconds\":", f);
	fprintf(f, "%.6f", t1 - t0);
	fputs("}", f);

	cert_free(&rpC);
	program_free(&rpP);
	if (ranrn) { cert_free(&rnC); program_free(&rnP); }
	free(col);
	free(cutS);
}

/* ---------------------------------------------------------------------- main */

static char *slurp(const char *path)
{
	FILE *f = path ? fopen(path, "rb") : stdin;
	size_t len = 0, cap = 1 << 16;
	char *buf;
	if (!f) { fprintf(stderr, "rn: cannot open %s\n", path); exit(2); }
	buf = xmalloc(cap);
	for (;;) {
		size_t k = fread(buf + len, 1, cap - len - 1, f);
		len += k;
		if (len + 1 < cap)
			break;
		cap *= 2;
		buf = xrealloc(buf, cap);
	}
	buf[len] = 0;
	if (path)
		fclose(f);
	return buf;
}

static const char usage[] =
"usage: rn [options] [file...]\n"
"\n"
"Decide whether each input graph is RP, SRN or not RN, and report a witness\n"
"point of the spanning tree polytope when there is one.  Input is JSON or\n"
"graph6 (auto-detected); with no file, stdin is read.\n"
"\n"
"  {\"name\": \"K_{2,3}\", \"vertices\": [\"a\",\"b\",...], \"edges\": [[\"a\",\"b\"], ...]}\n"
"  {\"name\": \"C6\", \"n\": 6, \"edges\": [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]]}\n"
"  {\"graph6\": \"D]o\"}        K_{2,3}; a bare \"D]o\" string works too\n"
"\n"
"A top-level array holds several graphs.\n"
"\n"
"options:\n"
"  --exact        settle the answer in rational arithmetic and report the\n"
"                 witness and the optimum exactly\n"
"  --jsonl        one result per line, with no enclosing array\n"
"  --no-witness   omit the witness point\n"
"  --no-echo      omit the vertex and edge lists from the output\n"
"  --no-graph6    omit the graph6 string from the output\n"
"  -h, --help     this message\n";

int main(int argc, char **argv)
{
	Opts o = { 1, 1, 1, 0, 0 };
	const char *files[256];
	int nf = 0, i, first = 1;
	FILE *f = stdout;

	glp_term_out(GLP_OFF);

	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
			{ fputs(usage, stdout); return 0; }
		else if (!strcmp(argv[i], "--exact")) o.exact = 1;
		else if (!strcmp(argv[i], "--jsonl")) o.jsonl = 1;
		else if (!strcmp(argv[i], "--no-witness")) o.witness = 0;
		else if (!strcmp(argv[i], "--no-echo")) o.echo = 0;
		else if (!strcmp(argv[i], "--no-graph6")) o.g6 = 0;
		else if (argv[i][0] == '-' && argv[i][1])
			{ fprintf(stderr, "rn: unknown option %s\n", argv[i]); return 2; }
		else if (nf < 256) files[nf++] = argv[i];
	}
	if (!nf)
		files[nf++] = NULL;

	if (!o.jsonl)
		fputs("[", f);

	for (i = 0; i < nf; i++) {
		char *buf = slurp(files[i]);
		const char *p = buf;

		while (*p && isspace((unsigned char)*p))
			p++;

		if (*p == '{' || *p == '[') {		/* JSON */
			JP j = { p, NULL };
			JV *root = jvalue(&j);
			int k, cnt;
			if (!root) {
				fprintf(stderr, "rn: %s: %s\n", files[i] ? files[i] : "(stdin)",
					j.err ? j.err : "parse error");
				free(buf);
				continue;
			}
			cnt = (root->t == JARR) ? root->nk : 1;
			for (k = 0; k < cnt; k++) {
				JV *e = (root->t == JARR) ? root->kid[k] : root;
				Graph g;
				const char *err;
				if (!graph_from_json(e, &g, &err)) {
					fprintf(stderr, "rn: graph %d: %s\n", k, err ? err : "bad entry");
					gfree(&g);
					continue;
				}
				if (!o.jsonl && !first) fputs(",\n", f);
				first = 0;
				run_one(f, &g, &o);
				if (o.jsonl) fputs("\n", f);
				gfree(&g);
			}
			jfree(root);
		} else {				/* graph6, one per line */
			char *line = buf, *nl;
			while (*line) {
				Graph g;
				nl = strchr(line, '\n');
				if (nl) *nl = 0;
				while (*line == ' ' || *line == '\t' || *line == '\r')
					line++;
				{
					char *e = line + strlen(line);
					while (e > line && (e[-1] == '\r' || e[-1] == ' '))
						*--e = 0;
				}
				if (*line) {
					memset(&g, 0, sizeof g);
					if (!g6_decode(line, &g)) {
						fprintf(stderr, "rn: bad graph6 line: %s\n", line);
					} else {
						normalise(&g);
						if (!o.jsonl && !first) fputs(",\n", f);
						first = 0;
						run_one(f, &g, &o);
						if (o.jsonl) fputs("\n", f);
					}
					gfree(&g);
				}
				if (!nl)
					break;
				line = nl + 1;
			}
		}
		free(buf);
	}

	if (!o.jsonl)
		fputs("]\n", f);
	return 0;
}
