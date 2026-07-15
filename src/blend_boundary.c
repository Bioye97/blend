/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

static void blend_boundary_free_points(double **points, int rows)
{
    int i;

    if (points == NULL) {
        return;
    }

    for (i = 0; i < rows; i++) {
        free(points[i]);
    }
    free(points);
}

static double **blend_boundary_copy_points(double **points, int rows)
{
    double **copy = NULL;
    int i;

    if (points == NULL || rows <= 0) {
        return NULL;
    }

    copy = (double **)malloc(sizeof(*copy) * (size_t)rows);
    if (copy == NULL) {
        return NULL;
    }

    for (i = 0; i < rows; i++) {
        copy[i] = (double *)malloc(sizeof(*copy[i]) * 2);
        if (copy[i] == NULL) {
            blend_boundary_free_points(copy, i);
            return NULL;
        }
        copy[i][0] = points[i][0];
        copy[i][1] = points[i][1];
    }

    return copy;
}

static void blend_boundary_restore_points(double **dst, double **src, int rows)
{
    int i;

    if (dst == NULL || src == NULL) {
        return;
    }

    for (i = 0; i < rows; i++) {
        dst[i][0] = src[i][0];
        dst[i][1] = src[i][1];
    }
}

void blend_permuted_vertex_free(permuted_vertex *vertex)
{
    if (vertex == NULL) {
        return;
    }

    blend_boundary_free_points(vertex->bottom_vertices, vertex->row_size_bv);
    blend_boundary_free_points(vertex->bottom_to_left_vertices, vertex->row_size_btlv);
    blend_boundary_free_points(vertex->left_vertices, vertex->row_size_lv);
    blend_boundary_free_points(vertex->left_to_top_vertices, vertex->row_size_lttv);
    blend_boundary_free_points(vertex->top_vertices, vertex->row_size_tv);
    blend_boundary_free_points(vertex->right_to_top_vertices, vertex->row_size_rttv);
    blend_boundary_free_points(vertex->right_vertices, vertex->row_size_rv);
    blend_boundary_free_points(vertex->bottom_to_right_vertices, vertex->row_size_btrv);
    memset(vertex, 0, sizeof(*vertex));
}

void blend_window_boundary_clear(window *data)
{
    if (data == NULL) {
        return;
    }

    blend_window_clear_polygon(data);
    free(data->nnx1);
    free(data->nnx2);
    free(data->nny1);
    free(data->nny2);
    data->nnx1 = NULL;
    data->nnx2 = NULL;
    data->nny1 = NULL;
    data->nny2 = NULL;
    data->minnx = 0;
    data->minny = 0;
}

int unique_vertices(double **array, int row_length) {
    /* Sort and update data->vertices and update data->row_size */

    /* Initialize local variables */
    int i, j, k, row;
    row = row_length;

    for (i = 0; i < row; i++) {
		for(j = i + 1; j < row;) {
    		if(array[i][0] == array[j][0] && array[i][1] == array[j][1]) {
    			for(k = j; k < row - 1; k++) {
    				array[k][0] = array[k+1][0];
                    array[k][1] = array[k+1][1];
				}
				row--;
			}
			else {
				j++;
			}
		}
	}

    /* Number of vertices error check */
    if (row < 3 ) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: the number of unique vertices is less than three; revise vertices\n");
		return FAIL;
	}

    return row;
}

int ascend_x(double **array, int row_length) {

    int i, j;
    double tmp_x, tmp_y;

    for (i = 0; i < row_length; i++) {
        for (j = i + 1; j < row_length; j++) {
            if (array[i][0] > array[j][0]) {
                tmp_x = array[i][0];
                tmp_y = array[i][1];
                array[i][0] = array[j][0];
                array[i][1] = array[j][1];
                array[j][0] = tmp_x;
                array[j][1] = tmp_y;
            }
        }
    }

    return SUCCESS;
}

int descend_x(double **array, int row_length) {

  int i, j;
    double tmp_x, tmp_y;

    for (i = 0; i < row_length; i++) {
        for (j = i + 1; j < row_length; j++) {
            if (array[i][0] < array[j][0]) {
                tmp_x = array[i][0];
                tmp_y = array[i][1];
                array[i][0] = array[j][0];
                array[i][1] = array[j][1];
                array[j][0] = tmp_x;
                array[j][1] = tmp_y;
            }
        }
    }
    return SUCCESS;
}

int ascend_y(double **array, int row_length) {

    int i, j;
    double tmp_x, tmp_y;

    for (i = 0; i < row_length; i++) {
        for (j = i + 1; j < row_length; j++) {
            if (array[i][1] > array[j][1]) {
                tmp_x = array[i][0];
                tmp_y = array[i][1];
                array[i][0] = array[j][0];
                array[i][1] = array[j][1];
                array[j][0] = tmp_x;
                array[j][1] = tmp_y;
            }
        }
    }

    return SUCCESS;
}

int descend_y(double **array, int row_length) {

  int i, j;
    double tmp_x, tmp_y;

    for (i = 0; i < row_length; i++) {
        for (j = i + 1; j < row_length; j++) {
            if (array[i][1] < array[j][1]) {
                tmp_x = array[i][0];
                tmp_y = array[i][1];
                array[i][0] = array[j][0];
                array[i][1] = array[j][1];
                array[j][0] = tmp_x;
                array[j][1] = tmp_y;
            }
        }
    }
    return SUCCESS;
}

int bounding_boxcheck(window *data) {

    /* Ensure that the range of the vertices is 
    [0 nx-1] and [0 ny-1] */

    /* Check that the polygon represented by the vertices 
    are bounded by the model dimension */

    /* Minimum bound check on x */
    if (ascend_x(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][0] < 0) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: minimum x-grid coordinate is less than 0\n");
        BLEND_Report(BLEND_MSG_ERROR, "boundary: a vertex is outside the model domain; revise vertices\n");
        return FAIL;
    }

    if (data->vertices[0][0] > 0) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: minimum x-grid coordinate is greater than 0\n");
        BLEND_Report(BLEND_MSG_ERROR, "boundary: polygon is not bounded by the model domain; revise vertices\n");
        return FAIL;
    }

    /* Maximum bound check on x */
     if (descend_x(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][0] < data->nx - 1) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: maximum x-grid coordinate is less than %d\n", data->nx-1);
        BLEND_Report(BLEND_MSG_ERROR, "boundary: polygon is not bounded by the model domain; revise vertices\n");
        return FAIL;
    }

    if (data->vertices[0][0] > data->nx - 1) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: maximum x-grid coordinate is greater than %d\n", data->nx-1);
        BLEND_Report(BLEND_MSG_ERROR, "boundary: a vertex is outside the model domain; revise vertices\n");
        return FAIL;
    }

    /* Minimum bound check on y */
    if (ascend_y(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][1] < 0) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: minimum y-grid coordinate is less than 0\n");
        BLEND_Report(BLEND_MSG_ERROR, "boundary: a vertex is outside the model domain; revise vertices\n");
        return FAIL;
    }

    if (data->vertices[0][1] > 0) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: minimum y-grid coordinate is greater than 0\n");
        BLEND_Report(BLEND_MSG_ERROR, "boundary: polygon is not bounded by the model domain; revise vertices\n");
        return FAIL;
    }

    /* Maximum bound check on y */
    if (descend_y(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][1] < data->ny - 1) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: maximum y-grid coordinate is less than %d\n", data->ny-1);
        BLEND_Report(BLEND_MSG_ERROR, "boundary: polygon is not bounded by the model domain; revise vertices\n");
        return FAIL;
    }

    if (data->vertices[0][1] > data->ny - 1) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: maximum y-grid coordinate is greater than %d\n", data->ny-1);
        BLEND_Report(BLEND_MSG_ERROR, "boundary: a vertex is outside the model domain; revise vertices\n");
        return FAIL;
    }

    return SUCCESS;
}

/* Assemble bottom vertices */
int assemble_bv(window *data, permuted_vertex *vertex) {

    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* Store row locations where y = 0 (minimum) */
    for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][1] == 0){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
    }

    /* Conditionals for assignment */
    if (count >= 2) {
        /* sort ascending in x and assign first and last rows */
        if (ascend_x(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* declare bottom vertices array */
	    vertex->bottom_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<2; m++){
            vertex->bottom_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* assign (min x, 0) (max x, 0) */
        vertex->bottom_vertices[0][0] = tmp_array[0][0];
        vertex->bottom_vertices[0][1] = tmp_array[0][1];
        vertex->bottom_vertices[1][0] = tmp_array[count-1][0];
        vertex->bottom_vertices[1][1] = tmp_array[count-1][1];

        /* assign number of bottom vertices */
        vertex->row_size_bv = 2;
    }
    /* Only one bottom vertex exists */
    else if (count == 1) {
       /* declare bottom vertices array */
	    vertex->bottom_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<1; m++){
            vertex->bottom_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }
        /* assign (min x, 0) (max x, 0) */
        vertex->bottom_vertices[0][0] = tmp_array[0][0];
        vertex->bottom_vertices[0][1] = tmp_array[0][1];
       /* assign number of bottom vertices */
       vertex->row_size_bv = 1;
    }
    else {
       BLEND_Report(BLEND_MSG_ERROR, "boundary: no bottom vertices found; revise vertices\n");
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble left vertices */
int assemble_lv(window *data, permuted_vertex *vertex) {
    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* store row locations where x = 0 (minimum) */
    for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] == 0){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
    }

    /* Conditionals for assignment */
    if (count >= 2) {
        /* sort ascending in y and assign first and last rows */
        if (ascend_y(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* declare left vertices array */
	    vertex->left_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<2; m++){
            vertex->left_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* Assign (0, min y) (0, max y) */
        vertex->left_vertices[0][0] = tmp_array[0][0];
        vertex->left_vertices[0][1] = tmp_array[0][1];
        vertex->left_vertices[1][0] = tmp_array[count-1][0];
        vertex->left_vertices[1][1] = tmp_array[count-1][1];

        /* assign number of left vertices */
        vertex->row_size_lv = 2;
    }
    /* Only one left vertex exists */
    else if (count == 1) {
       /* declare left vertices array */
	    vertex->left_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<1; m++){
            vertex->left_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }
        /* Assign (0, min y) (0, max y) */
        vertex->left_vertices[0][0] = tmp_array[0][0];
        vertex->left_vertices[0][1] = tmp_array[0][1];
       /* assign number of left vertices */
       vertex->row_size_lv = 1;
    }
    else {
       BLEND_Report(BLEND_MSG_ERROR, "boundary: no left vertices found; revise vertices\n");
    }

    free(tmp_array);
    return SUCCESS;

}

/* Assemble top vertices */
int assemble_tv(window *data, permuted_vertex *vertex) {
    
    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* store row locations where y = ny - 1 (maximum) */
    for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][1] == data->ny-1) {
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
    }

    /* Conditionals for assignment */
    if (count >= 2) {
        /* sort ascending in x and assign first and last rows */
        if (ascend_x(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* declare bottom vertices array */
	    vertex->top_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<2; m++){
            vertex->top_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* Assign (min x, ny-1) (max x, ny-1) */
        vertex->top_vertices[0][0] = tmp_array[0][0];
        vertex->top_vertices[0][1] = tmp_array[0][1];
        vertex->top_vertices[1][0] = tmp_array[count-1][0];
        vertex->top_vertices[1][1] = tmp_array[count-1][1];

        /* assign number of bottom vertices */
        vertex->row_size_tv = 2;
    }
    /* Only one bottom vertex exists */
    else if (count == 1) {
       /* declare bottom vertices array */
	    vertex->top_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<1; m++){
            vertex->top_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }
        /* Assign (min x, ny-1) (max x, ny-1) */
        vertex->top_vertices[0][0] = tmp_array[0][0];
        vertex->top_vertices[0][1] = tmp_array[0][1];
       /* assign number of bottom vertices */
       vertex->row_size_tv = 1;
    }
    else {
       BLEND_Report(BLEND_MSG_ERROR, "boundary: no top vertices found; revise vertices\n");
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble right vertices */
int assemble_rv(window *data, permuted_vertex *vertex) {
    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* store row locations where x = nx-1 (maximum) */
    for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] == data->nx-1){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
    }

    /* Conditionals for assignment */
    if (count >= 2) {
        /* sort ascending in y and assign first and last rows */
        if (ascend_y(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* declare right vertices array */
	    vertex->right_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<2; m++){
            vertex->right_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* Assign (nx-1, min y) (nx-1, max y) */
        vertex->right_vertices[0][0] = tmp_array[0][0];
        vertex->right_vertices[0][1] = tmp_array[0][1];
        vertex->right_vertices[1][0] = tmp_array[count-1][0];
        vertex->right_vertices[1][1] = tmp_array[count-1][1];

        /* assign number of left vertices */
        vertex->row_size_rv = 2;
    }
    /* Only one left vertex exists */
    else if (count == 1) {
       /* declare left vertices array */
	    vertex->right_vertices = (double **)malloc(sizeof(double *)* 2);
        for (m=0; m<1; m++){
            vertex->right_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }
        /* Assign (nx-1, min y) (nx-1, max y) */
        vertex->right_vertices[0][0] = tmp_array[0][0];
        vertex->right_vertices[0][1] = tmp_array[0][1];
       /* assign number of left vertices */
       vertex->row_size_rv = 1;
    }
    else {
       BLEND_Report(BLEND_MSG_ERROR, "boundary: no right vertices found; revise vertices\n");
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble bottom-to-left-vertices */
int assemble_btlv(window *data, permuted_vertex *vertex) {

    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* store row locations where x and y are between bottom_vertices[1]
    and left_vertices[1] (maximum) */
    for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] < vertex->bottom_vertices[0][0] && data->vertices[i][0] > vertex->left_vertices[0][0] 
        &&  data->vertices[i][1] > vertex->bottom_vertices[0][1] && data->vertices[i][1] < vertex->left_vertices[0][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
    }

    /* Conditionals for assignment */
    if (count >= 1) {
        /* sort ascending in y */
        if (ascend_y(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* Check for irregular vertices, i.e., y should be strictly increasing and
        x should be strictly decrasing in this range. */
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][1] > tmp_array[m+1][1]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from bottom edge to left edge; y coordinates should be strictly increasing\n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] < tmp_array[m+1][0]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from bottom edge to left edge; x coordinates should be strictly decreasing\n");
                return FAIL;
            }
        }

        /* declare bottom-to-left vertices array */
	    vertex->bottom_to_left_vertices = (double **)malloc(sizeof(double *)* count);
        for (m=0; m<count; m++){
            vertex->bottom_to_left_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* Assignment */
        for (m=0; m<count; m++){
            vertex->bottom_to_left_vertices[m][0] = tmp_array[m][0];
            vertex->bottom_to_left_vertices[m][1] = tmp_array[m][1];
        }
        
        /* assign number of left vertices */
        vertex->row_size_btlv = count;
    }
    else {
       /* No bottom-to-left vertices */ 
       vertex->row_size_btlv = 0;
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble left-to-top vertices */
int assemble_lttv(window *data, permuted_vertex *vertex){

    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* store row locations where x and y are between left_vertices[1 or 2]
    and top_vertices[1] (maximum) */
    if (vertex->row_size_lv < 2) {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] > vertex->left_vertices[0][0] && data->vertices[i][0] < vertex->top_vertices[0][0] 
        &&  data->vertices[i][1] > vertex->left_vertices[0][1] && data->vertices[i][1] < vertex->top_vertices[0][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }
    }

    else {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] > vertex->left_vertices[1][0] && data->vertices[i][0] < vertex->top_vertices[0][0] 
        &&  data->vertices[i][1] > vertex->left_vertices[1][1] && data->vertices[i][1] < vertex->top_vertices[0][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }
    }
    
    /* Conditionals for assignment */
    if (count >= 1) {
        /* sort ascending in y */
        if (ascend_y(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* Check for irregular vertices, i.e., y should be strictly increasing and
        x should be strictly increasing in this range. */
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][1] > tmp_array[m+1][1]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from left edge to top edge; y coordinates should be strictly increasing\n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] > tmp_array[m+1][0]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from left edge to top edge; x coordinates should be strictly increasing\n");
                return FAIL;
            }
        }

        /* declare left-to-top vertices array */
	    vertex->left_to_top_vertices = (double **)malloc(sizeof(double *)* count);
        for (m=0; m<count; m++){
            vertex->left_to_top_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* Assignment */
        for (m=0; m<count; m++){
            vertex->left_to_top_vertices[m][0] = tmp_array[m][0];
            vertex->left_to_top_vertices[m][1] = tmp_array[m][1];
        }
        
        /* assign number of left-to-top vertices */
        vertex->row_size_lttv = count;
    }
    else {
       /* No left-to-top vertices */ 
       vertex->row_size_lttv = 0;
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble right-to-top vertices */
int assemble_rttv(window *data, permuted_vertex *vertex) {

    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* store row locations where x and y are between right_vertices[1 or 2]
    and top_vertices[1 or 2] */
    if (vertex->row_size_rv < 2 && vertex->row_size_tv < 2) {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] < vertex->right_vertices[0][0] && data->vertices[i][0] > vertex->top_vertices[0][0] 
        &&  data->vertices[i][1] > vertex->right_vertices[0][1] && data->vertices[i][1] < vertex->top_vertices[0][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }
    }

    else if (vertex->row_size_rv == 2 && vertex->row_size_tv < 2) {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] < vertex->right_vertices[1][0] && data->vertices[i][0] > vertex->top_vertices[0][0] 
        &&  data->vertices[i][1] > vertex->right_vertices[1][1] && data->vertices[i][1] < vertex->top_vertices[0][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }
    }

    else if (vertex->row_size_rv < 2 && vertex->row_size_tv == 2) {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] < vertex->right_vertices[0][0] && data->vertices[i][0] > vertex->top_vertices[1][0] 
        &&  data->vertices[i][1] > vertex->right_vertices[0][1] && data->vertices[i][1] < vertex->top_vertices[1][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }
    }

    /* two vertices exist here */
    else {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] < vertex->right_vertices[1][0] && data->vertices[i][0] > vertex->top_vertices[1][0] 
        &&  data->vertices[i][1] > vertex->right_vertices[1][1] && data->vertices[i][1] < vertex->top_vertices[1][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }

    }
    
    /* Conditionals for assignment */
    if (count >= 1) {
        /* sort ascending in y */
        if (ascend_y(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* Check for irregular vertices, i.e., y should be strictly increasing and
        x should be strictly decreasing in this range. */
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][1] > tmp_array[m+1][1]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from right edge to top edge; y coordinates should be strictly increasing\n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] < tmp_array[m+1][0]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from right edge to top edge; x coordinates should be strictly decreasing\n");
                return FAIL;
            }
        }

        /* declare right-to-top vertices array */
	    vertex->right_to_top_vertices = (double **)malloc(sizeof(double *)* count);
        for (m=0; m<count; m++){
            vertex->right_to_top_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* Assignment */
        for (m=0; m<count; m++){
            vertex->right_to_top_vertices[m][0] = tmp_array[m][0];
            vertex->right_to_top_vertices[m][1] = tmp_array[m][1];
        }
        
        /* assign number of right-to-top vertices */
        vertex->row_size_rttv = count;
    }
    else {
       /* No right-to-top vertices */ 
       vertex->row_size_rttv = 0;
    }

    free(tmp_array);
    return SUCCESS;

}

/* Assemble bottom-to-right vertices */
int assemble_btrv(window *data, permuted_vertex *vertex) {
        
    int i, m, j=0, count=0;
    const int length = data->row_size;
    double **tmp_array;
	tmp_array = (double **)malloc(sizeof(double *)* length);
    for (m=0; m<length; m++){
        tmp_array[m] = (double *)malloc(sizeof(double)* 2); 
    }

    /* store row locations where x and y are between bottom_vertices[1 or 2]
    and right_vertices[1] (maximum) */

    if (vertex->row_size_bv < 2) {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] > vertex->bottom_vertices[0][0] && data->vertices[i][0] < vertex->right_vertices[0][0] 
        &&  data->vertices[i][1] > vertex->bottom_vertices[0][1] && data->vertices[i][1] < vertex->right_vertices[0][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }
    }
    else {

        for (i=0; i<data->row_size; i++) {
        if (data->vertices[i][0] > vertex->bottom_vertices[1][0] && data->vertices[i][0] < vertex->right_vertices[0][0] 
        &&  data->vertices[i][1] > vertex->bottom_vertices[1][1] && data->vertices[i][1] < vertex->right_vertices[0][1]){
            /* store rows */
            tmp_array[j][0] = data->vertices[i][0];
            tmp_array[j][1] = data->vertices[i][1];
            /* store size */
            count++;
            /* update index */
            j++;
        }
        }
    }

    /* Conditionals for assignment */
    if (count >= 1) {
        /* sort ascending in y */
        if (ascend_y(tmp_array, count) != SUCCESS) {
            return FAIL;
        }

        /* Check for irregular vertices, i.e., y should be strictly increasing and
        x should be strictly increasing in this range. */
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][1] > tmp_array[m+1][1]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from bottom edge to right edge; y coordinates should be strictly increasing\n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] > tmp_array[m+1][0]) {
                BLEND_Report(BLEND_MSG_ERROR, "boundary: irregular vertices from bottom edge to right edge; x coordinates should be strictly increasing\n");
                return FAIL;
            }
        }

        /* declare bottom-to-right vertices array */
	    vertex->bottom_to_right_vertices = (double **)malloc(sizeof(double *)* count);
        for (m=0; m<count; m++){
            vertex->bottom_to_right_vertices[m] = (double *)malloc(sizeof(double)* 2); 
        }

        /* Assignment */
        for (m=0; m<count; m++){
            vertex->bottom_to_right_vertices[m][0] = tmp_array[m][0];
            vertex->bottom_to_right_vertices[m][1] = tmp_array[m][1];
        }
        
        /* assign number of bottom-to-right vertices */
        vertex->row_size_btrv = count;
    }
    else {
       /* No bottom-to-right vertices */ 
       vertex->row_size_btrv = 0;
    }

    free(tmp_array);
    return SUCCESS;
}

typedef struct boundary_path_candidate {
    int valid;
    int count;
    int *indices;
    double *points;
} boundary_path_candidate;

static void boundary_path_candidate_free(boundary_path_candidate *candidate)
{
    if (candidate == NULL) {
        return;
    }

    free(candidate->indices);
    free(candidate->points);
    candidate->indices = NULL;
    candidate->points = NULL;
    candidate->count = 0;
    candidate->valid = 0;
}

static int boundary_path_candidate_alloc(boundary_path_candidate *candidate, int capacity)
{
    if (candidate == NULL || capacity < 0) {
        return FAIL;
    }

    memset(candidate, 0, sizeof(*candidate));
    candidate->indices = (int *)calloc((size_t)(capacity > 0 ? capacity : 1), sizeof(*candidate->indices));
    candidate->points = (double *)calloc((size_t)(capacity > 0 ? capacity : 1) * 2, sizeof(*candidate->points));
    if (candidate->indices == NULL || candidate->points == NULL) {
        boundary_path_candidate_free(candidate);
        BLEND_Report(BLEND_MSG_ERROR, "boundary: could not allocate traversal path candidate\n");
        return FAIL;
    }

    return SUCCESS;
}

static double boundary_coordinate_tolerance(const window *data)
{
    double scale = 1.0;

    if (data != NULL) {
        if (data->nx > 0) scale = fmax(scale, (double)data->nx);
        if (data->ny > 0) scale = fmax(scale, (double)data->ny);
    }

    return 1.0e-10 * scale;
}

static int boundary_points_equal(const double *a, const double *b, double tolerance)
{
    return fabs(a[0] - b[0]) <= tolerance && fabs(a[1] - b[1]) <= tolerance;
}

static int boundary_value_between_open(double value, double a, double b, double tolerance)
{
    double min_value = fmin(a, b);
    double max_value = fmax(a, b);

    return value > min_value + tolerance && value < max_value - tolerance;
}

static int boundary_point_on_bbox_side(const window *data, const double *point, double tolerance)
{
    return fabs(point[0]) <= tolerance ||
           fabs(point[0] - (double)(data->nx - 1)) <= tolerance ||
           fabs(point[1]) <= tolerance ||
           fabs(point[1] - (double)(data->ny - 1)) <= tolerance;
}

static int boundary_find_vertex_index(const window *data, const double *point, double tolerance, int *index_out)
{
    int i;

    if (data == NULL || point == NULL || index_out == NULL) {
        return FAIL;
    }

    for (i = 0; i < data->row_size; i++) {
        if (fabs(data->vertices[i][0] - point[0]) <= tolerance &&
            fabs(data->vertices[i][1] - point[1]) <= tolerance) {
            *index_out = i;
            return SUCCESS;
        }
    }

    return FAIL;
}

static int boundary_path_is_strictly_monotone(const double *start, const boundary_path_candidate *candidate,
                                              const double *end, int x_direction, int y_direction,
                                              double tolerance)
{
    double prev_x = start[0];
    double prev_y = start[1];
    int i;

    for (i = 0; i < candidate->count; i++) {
        double x = candidate->points[2 * i];
        double y = candidate->points[2 * i + 1];
        double dx = x - prev_x;
        double dy = y - prev_y;

        if ((double)x_direction * dx <= tolerance || (double)y_direction * dy <= tolerance) {
            return 0;
        }
        prev_x = x;
        prev_y = y;
    }

    if ((double)x_direction * (end[0] - prev_x) <= tolerance ||
        (double)y_direction * (end[1] - prev_y) <= tolerance) {
        return 0;
    }

    return 1;
}

static int boundary_collect_path_candidate(const window *data, int start_index, int end_index, int step,
                                           const double *start, const double *end,
                                           int x_direction, int y_direction, double tolerance,
                                           boundary_path_candidate *candidate)
{
    int index;
    int guard = 0;

    if (data == NULL || start == NULL || end == NULL || candidate == NULL ||
        (step != 1 && step != -1) || data->row_size <= 0) {
        return FAIL;
    }

    candidate->valid = 0;
    candidate->count = 0;

    if (start_index == end_index) {
        candidate->valid = 1;
        return SUCCESS;
    }

    index = start_index;
    while (index != end_index) {
        double *point;

        if (step > 0) {
            index = (index + 1) % data->row_size;
        }
        else {
            index = (index == 0) ? data->row_size - 1 : index - 1;
        }
        guard++;
        if (guard > data->row_size) {
            return FAIL;
        }
        if (index == end_index) {
            break;
        }

        point = data->vertices[index];
        if (boundary_point_on_bbox_side(data, point, tolerance)) {
            continue;
        }
        if (!boundary_value_between_open(point[0], start[0], end[0], tolerance) ||
            !boundary_value_between_open(point[1], start[1], end[1], tolerance)) {
            return SUCCESS;
        }

        candidate->indices[candidate->count] = index;
        candidate->points[2 * candidate->count] = point[0];
        candidate->points[2 * candidate->count + 1] = point[1];
        candidate->count++;
    }

    if (candidate->count == 0) {
        if (boundary_points_equal(start, end, tolerance) ||
            ((double)x_direction * (end[0] - start[0]) > tolerance &&
             (double)y_direction * (end[1] - start[1]) > tolerance)) {
            candidate->valid = 1;
        }
        return SUCCESS;
    }

    candidate->valid = boundary_path_is_strictly_monotone(start, candidate, end,
                                                          x_direction, y_direction, tolerance);
    return SUCCESS;
}

static int boundary_nonbbox_vertex_count(const window *data, double tolerance)
{
    int i;
    int count = 0;

    for (i = 0; i < data->row_size; i++) {
        if (!boundary_point_on_bbox_side(data, data->vertices[i], tolerance)) {
            count++;
        }
    }

    return count;
}

static int boundary_choose_path_candidates(const window *data, double tolerance,
                                           boundary_path_candidate candidates[4][2],
                                           int selected[4])
{
    int a, b, c, d;
    int nonbbox_count = boundary_nonbbox_vertex_count(data, tolerance);

    for (a = 0; a < 2; a++) {
        for (b = 0; b < 2; b++) {
            for (c = 0; c < 2; c++) {
                for (d = 0; d < 2; d++) {
                    int choices[4];
                    unsigned char *assigned = NULL;
                    int assigned_count = 0;
                    int sector;
                    int ok = 1;

                    choices[0] = a;
                    choices[1] = b;
                    choices[2] = c;
                    choices[3] = d;

                    for (sector = 0; sector < 4; sector++) {
                        if (!candidates[sector][choices[sector]].valid) {
                            ok = 0;
                            break;
                        }
                    }
                    if (!ok) {
                        continue;
                    }

                    assigned = (unsigned char *)calloc((size_t)(data->row_size > 0 ? data->row_size : 1),
                                                       sizeof(*assigned));
                    if (assigned == NULL) {
                        BLEND_Report(BLEND_MSG_ERROR, "boundary: could not allocate traversal assignments\n");
                        return FAIL;
                    }

                    for (sector = 0; sector < 4 && ok; sector++) {
                        boundary_path_candidate *candidate = &candidates[sector][choices[sector]];
                        int i;

                        for (i = 0; i < candidate->count; i++) {
                            int index = candidate->indices[i];

                            if (assigned[index]) {
                                ok = 0;
                                break;
                            }
                            assigned[index] = 1;
                            assigned_count++;
                        }
                    }

                    if (ok && assigned_count == nonbbox_count) {
                        selected[0] = choices[0];
                        selected[1] = choices[1];
                        selected[2] = choices[2];
                        selected[3] = choices[3];
                        free(assigned);
                        return SUCCESS;
                    }

                    free(assigned);
                }
            }
        }
    }

    return FAIL;
}

static int boundary_assign_path_vertices(double ***target, int *row_size,
                                         const boundary_path_candidate *candidate)
{
    double **points = NULL;
    int i;

    if (target == NULL || row_size == NULL || candidate == NULL) {
        return FAIL;
    }

    *target = NULL;
    *row_size = 0;

    if (candidate->count == 0) {
        return SUCCESS;
    }

    points = (double **)malloc(sizeof(*points) * (size_t)candidate->count);
    if (points == NULL) {
        return FAIL;
    }
    for (i = 0; i < candidate->count; i++) {
        points[i] = (double *)malloc(sizeof(*points[i]) * 2);
        if (points[i] == NULL) {
            blend_boundary_free_points(points, i);
            return FAIL;
        }
        points[i][0] = candidate->points[2 * i];
        points[i][1] = candidate->points[2 * i + 1];
    }

    if (ascend_y(points, candidate->count) != SUCCESS) {
        blend_boundary_free_points(points, candidate->count);
        return FAIL;
    }

    *target = points;
    *row_size = candidate->count;
    return SUCCESS;
}

static int assemble_path_vertices(window *data, permuted_vertex *vertex)
{
    boundary_path_candidate candidates[4][2];
    double tolerance = boundary_coordinate_tolerance(data);
    double *bottom_left = vertex->bottom_vertices[0];
    double *bottom_right = vertex->bottom_vertices[vertex->row_size_bv > 1 ? 1 : 0];
    double *left_bottom = vertex->left_vertices[0];
    double *left_top = vertex->left_vertices[vertex->row_size_lv > 1 ? 1 : 0];
    double *top_left = vertex->top_vertices[0];
    double *top_right = vertex->top_vertices[vertex->row_size_tv > 1 ? 1 : 0];
    double *right_bottom = vertex->right_vertices[0];
    double *right_top = vertex->right_vertices[vertex->row_size_rv > 1 ? 1 : 0];
    int bottom_left_index = 0, bottom_right_index = 0;
    int left_bottom_index = 0, left_top_index = 0;
    int top_left_index = 0, top_right_index = 0;
    int right_bottom_index = 0, right_top_index = 0;
    int selected[4] = {-1, -1, -1, -1};
    int i;
    int initialized = 0;
    int status = FAIL;

    memset(candidates, 0, sizeof(candidates));
    for (i = 0; i < 8; i++) {
        if (boundary_path_candidate_alloc(&candidates[i / 2][i % 2], data->row_size) != SUCCESS) {
            goto cleanup;
        }
    }
    initialized = 1;

    if (boundary_find_vertex_index(data, bottom_left, tolerance, &bottom_left_index) != SUCCESS ||
        boundary_find_vertex_index(data, bottom_right, tolerance, &bottom_right_index) != SUCCESS ||
        boundary_find_vertex_index(data, left_bottom, tolerance, &left_bottom_index) != SUCCESS ||
        boundary_find_vertex_index(data, left_top, tolerance, &left_top_index) != SUCCESS ||
        boundary_find_vertex_index(data, top_left, tolerance, &top_left_index) != SUCCESS ||
        boundary_find_vertex_index(data, top_right, tolerance, &top_right_index) != SUCCESS ||
        boundary_find_vertex_index(data, right_bottom, tolerance, &right_bottom_index) != SUCCESS ||
        boundary_find_vertex_index(data, right_top, tolerance, &right_top_index) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "boundary: could not locate bounding-box anchor vertices\n");
        goto cleanup;
    }

    if (boundary_collect_path_candidate(data, bottom_left_index, left_bottom_index, 1,
                                        bottom_left, left_bottom, -1, 1, tolerance,
                                        &candidates[0][0]) != SUCCESS ||
        boundary_collect_path_candidate(data, bottom_left_index, left_bottom_index, -1,
                                        bottom_left, left_bottom, -1, 1, tolerance,
                                        &candidates[0][1]) != SUCCESS ||
        boundary_collect_path_candidate(data, left_top_index, top_left_index, 1,
                                        left_top, top_left, 1, 1, tolerance,
                                        &candidates[1][0]) != SUCCESS ||
        boundary_collect_path_candidate(data, left_top_index, top_left_index, -1,
                                        left_top, top_left, 1, 1, tolerance,
                                        &candidates[1][1]) != SUCCESS ||
        boundary_collect_path_candidate(data, right_top_index, top_right_index, 1,
                                        right_top, top_right, -1, 1, tolerance,
                                        &candidates[2][0]) != SUCCESS ||
        boundary_collect_path_candidate(data, right_top_index, top_right_index, -1,
                                        right_top, top_right, -1, 1, tolerance,
                                        &candidates[2][1]) != SUCCESS ||
        boundary_collect_path_candidate(data, bottom_right_index, right_bottom_index, 1,
                                        bottom_right, right_bottom, 1, 1, tolerance,
                                        &candidates[3][0]) != SUCCESS ||
        boundary_collect_path_candidate(data, bottom_right_index, right_bottom_index, -1,
                                        bottom_right, right_bottom, 1, 1, tolerance,
                                        &candidates[3][1]) != SUCCESS) {
        goto cleanup;
    }

    if (boundary_choose_path_candidates(data, tolerance, candidates, selected) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR,
                     "boundary: polygon vertices cannot be assigned to strict xy-monotone traversal paths; revise vertices\n");
        goto cleanup;
    }

    if (boundary_assign_path_vertices(&vertex->bottom_to_left_vertices, &vertex->row_size_btlv,
                                      &candidates[0][selected[0]]) != SUCCESS ||
        boundary_assign_path_vertices(&vertex->left_to_top_vertices, &vertex->row_size_lttv,
                                      &candidates[1][selected[1]]) != SUCCESS ||
        boundary_assign_path_vertices(&vertex->right_to_top_vertices, &vertex->row_size_rttv,
                                      &candidates[2][selected[2]]) != SUCCESS ||
        boundary_assign_path_vertices(&vertex->bottom_to_right_vertices, &vertex->row_size_btrv,
                                      &candidates[3][selected[3]]) != SUCCESS) {
        goto cleanup;
    }

    status = SUCCESS;

cleanup:
    if (status != SUCCESS) {
        blend_boundary_free_points(vertex->bottom_to_left_vertices, vertex->row_size_btlv);
        blend_boundary_free_points(vertex->left_to_top_vertices, vertex->row_size_lttv);
        blend_boundary_free_points(vertex->right_to_top_vertices, vertex->row_size_rttv);
        blend_boundary_free_points(vertex->bottom_to_right_vertices, vertex->row_size_btrv);
        vertex->bottom_to_left_vertices = NULL;
        vertex->left_to_top_vertices = NULL;
        vertex->right_to_top_vertices = NULL;
        vertex->bottom_to_right_vertices = NULL;
        vertex->row_size_btlv = 0;
        vertex->row_size_lttv = 0;
        vertex->row_size_rttv = 0;
        vertex->row_size_btrv = 0;
    }
    if (initialized) {
        for (i = 0; i < 8; i++) {
            boundary_path_candidate_free(&candidates[i / 2][i % 2]);
        }
    }
    return status;
}

/* Hanging sweep */
int hanging_sweep(int j, double *vector, double **prev_array, int prev_row, double **sub_array, int sub_row,
                  double **hang_array, int hang_len,  assembly_type step_type) {

    int index, i=0;
    double delta;

    /* Hanging array has elements */
    if (hang_len > 0) {

        /* Begin switch */
        switch(step_type) {

             /* ROW_STEP case */
             case ROW_STEP:

                /* Begin loop */
                for (index=j; index<sub_array[sub_row][1]; index++) {

                    /* index = j */
                    if (index == j){

                        /* compute delta */
                        delta = (prev_array[prev_row][0] - hang_array[i][0])/
                                (prev_array[prev_row][1] - hang_array[i][1]);

                        /* update vector with first member from hanging vertices */
                        vector[index] = vector[index-1] + delta;
                    }

                    /* index = hang_array[i][1] && hang_len = 1: 
                    one hanging vertex */
                    else if (index == hang_array[i][1] && hang_len == 1) {

                        /* Assign current value */
                        vector[index] = hang_array[i][0];

                        /* compute delta using the subsequent edge vertex */
                        delta = (hang_array[i][0] - sub_array[sub_row][0])/
                                (hang_array[i][1] - sub_array[sub_row][1]);

                    }

                    /* index = hang_array[i][1] && hang_len > 1 && i != hang_len: 
                    update between hanging vertices */
                    else if (index == hang_array[i][1] && hang_len > 1 && i != hang_len - 1) {

                        /* Assign current value */
                        vector[index] = hang_array[i][0];

                        /* compute delta using the next member in the hanging array */
                        delta = (hang_array[i][0] - hang_array[i+1][0])/
                                (hang_array[i][1] - hang_array[i+1][1]);

                        /* update hanging array index */
                        i++;

                    }

                    /* index = hang_array[hang_len-1][1]: last hanging vertex */
                    else if (index == hang_array[hang_len-1][1]) {

                        /* Assign value to vector */
                        vector[index] = hang_array[hang_len-1][0];

                        /* compute delta using the subsequent edge vertex */
                        delta = (hang_array[i][0] - sub_array[sub_row][0])/
                                (hang_array[i][1] - sub_array[sub_row][1]);
                    }

                    /* spatial marching with delta */
                    else {
                        vector[index] = vector[index-1] + delta;
                    }
                }
                index--; /* update index for next vertex condition */
                break;


             /* COLUMN_STEP case */
             case COLUMN_STEP:
                /* Begin loop */
                for (index=j; index<sub_array[sub_row][0]; index++) {

                    /* index = j */
                    if (index == j){

                        /* compute delta */
                        delta = (prev_array[prev_row][1] - hang_array[i][1])/
                                (prev_array[prev_row][0] - hang_array[i][0]);

                        /* update vector with first member from hanging vertices */
                        vector[index] = vector[index-1] + delta;
                    }

                    /* index = hang_array[i][1] && hang_len = 1: 
                    one hanging vertex */
                    else if (index == hang_array[i][0] && hang_len == 1) {

                        /* Assign current value */
                        vector[index] = hang_array[i][1];

                        /* compute delta using the subsequent edge vertex */
                        delta = (hang_array[i][1] - sub_array[sub_row][1])/
                                (hang_array[i][0] - sub_array[sub_row][0]);

                    }

                    /* index = hang_array[i][1] && hang_len > 1 && i != hang_len: 
                    update between hanging vertices */
                    else if (index == hang_array[i][0] && hang_len > 1 && i != hang_len - 1) {

                        /* Assign current value */
                        vector[index] = hang_array[i][1];

                        /* compute delta using the next member in the hanging array */
                        delta = (hang_array[i][1] - hang_array[i+1][1])/
                                (hang_array[i][0] - hang_array[i+1][0]);

                        /* update hanging array index */
                        i++;

                    }

                    /* index = hang_array[hang_len-1][1]: last hanging vertex */
                    else if (index == hang_array[hang_len-1][0]) {

                        /* Assign value to vector */
                        vector[index] = hang_array[hang_len-1][1];

                        /* compute delta using the subsequent edge vertex */
                        delta = (hang_array[i][1] - sub_array[sub_row][1])/
                                (hang_array[i][0] - sub_array[sub_row][0]);
                    }

                    /* spatial marching with delta */
                    else {
                        vector[index] = vector[index-1] + delta;
                    }
                }
                index--; /* update index for next vertex condition */
                break;
        } /* end of switch block */
    } 
    /* Hanging array has no elements, so index is the same */
    else {

        index = j;
    }

    return index;
}


/* Assemble minimum row-stepping vectors: nnx1 */
int minimum_rowstep(window *data, permuted_vertex *vertex) {

    /* bottom -> left -> top */
    int j=0; 
    double delta_x1;
    double *nnx1;
    const int length = data->ny;
    /* Initialize the nnx1 vector: dimension is ny */
	nnx1 = (double *)malloc(sizeof(double) * length);
    data->nnx1 = (int *)malloc(sizeof(int) * length);
    /* stepping mode */
    data->amode = ROW_STEP;

    /* first make sure the vertices are sorted correctly */
    if (ascend_x(vertex->bottom_vertices, vertex->row_size_bv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->bottom_to_left_vertices, vertex->row_size_btlv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->left_vertices, vertex->row_size_lv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->left_to_top_vertices, vertex->row_size_lttv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_x(vertex->top_vertices, vertex->row_size_tv) != SUCCESS) {
            return FAIL;
    }

    /* Begin */
    while (j < data->ny) {

        /* condition at first bottom edge */
        if (j == vertex->bottom_vertices[0][1]) {

            /* First assign first element to nnx1, i.e., x component of bottom_vertices */
            nnx1[j] = vertex->bottom_vertices[0][0];

            /* ~ Directions for stepping to the next edge ~ */

            /* Corner vertex exists : bottom and left corner */
            if (vertex->bottom_vertices[0][1] == vertex->left_vertices[0][1]) {

                /* Second left vertex does not exist */
                if (vertex->row_size_lv < 2) {

                    /* Left-to-top vertices exist */
                    if (vertex->row_size_lttv > 0) {

                        /* Call hanging_sweep to proceed */
                        j = hanging_sweep(j+1, nnx1, vertex->left_vertices, 0, vertex->top_vertices, 0,
                                          vertex->left_to_top_vertices, vertex->row_size_lttv, data->amode);
                    }

                    /* No left-to-top vertices */
                    else {
                        /* compute delta */
                        delta_x1 = (vertex->left_vertices[0][0] - vertex->top_vertices[0][0])/
                                   (vertex->left_vertices[0][1] - vertex->top_vertices[0][1]);
                    }
                }

                /* Second left vertex exits */
                else {
                    /* compute delta */
                    delta_x1 = (vertex->left_vertices[0][0] - vertex->left_vertices[1][0])/
                               (vertex->left_vertices[0][1] - vertex->left_vertices[1][1]);
                }

            }

            /* No corner vertex but bottom-to-left vertices exist */
            else if (vertex->row_size_btlv > 0) {

                /* Call hanging_sweep to proceed */
                j = hanging_sweep(j+1, nnx1, vertex->bottom_vertices, 0, vertex->left_vertices, 0,
                                  vertex->bottom_to_left_vertices, vertex->row_size_btlv, data->amode);
                
            }

            /* No corner vertex and no bottom-to-left vertices */
            else {
                /* compute delta */
                delta_x1 = (vertex->bottom_vertices[0][0] - vertex->left_vertices[0][0])/
                           (vertex->bottom_vertices[0][1] - vertex->left_vertices[0][1]);
            }

        }


        /* condition at first left edge */
        else if (j == vertex->left_vertices[0][1]) {

            /* Assign nnx1 value */
            nnx1[j] = vertex->left_vertices[0][0];

            /* ~ Directions for stepping to the next edge ~ */

            /* Second left vertex does not exist */
            if (vertex->row_size_lv < 2) {

                /* left-to-top vertices exist */
                if (vertex->row_size_lttv > 0) {
                    /* Call haanging_sweep to proceed */
                    j = hanging_sweep(j+1, nnx1, vertex->left_vertices, 0, vertex->top_vertices, 0,
                                          vertex->left_to_top_vertices, vertex->row_size_lttv, data->amode);

                }

                /* No left-to_top vertices */
                else {
                    /* compute delta */
                    delta_x1 = (vertex->left_vertices[0][0] - vertex->top_vertices[0][0])/
                               (vertex->left_vertices[0][1] - vertex->top_vertices[0][1]);
                }

            }

            /* Second left vertex exists */
            else {
                /* compute delta */
                delta_x1 = (vertex->left_vertices[0][0] - vertex->left_vertices[1][0])/
                           (vertex->left_vertices[0][1] - vertex->left_vertices[1][1]);
            }
        }

        /* condition at second left edge */
        else if (vertex->row_size_lv == 2 && j == vertex->left_vertices[1][1]) {

            /* Corner vertex exists */
            if (vertex->left_vertices[1][1] == vertex->top_vertices[0][1]) {

                /* Just assign nnx1 and we are done */
                nnx1[j] = vertex->left_vertices[1][0];
                break;
            }

            /* No corner vertex */
            else {

                /* Assign nnx1 */
                nnx1[j] = vertex->left_vertices[1][0];

                /* Left-to-top vertices exist */
                if (vertex->row_size_lttv > 0) {

                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nnx1, vertex->left_vertices, 1, vertex->top_vertices, 0,
                                      vertex->left_to_top_vertices, vertex->row_size_lttv, data->amode);
                    

                }

                /* No left-to-top vertices */
                else {
                    /* compute delta */
                    delta_x1 = (vertex->left_vertices[1][0] - vertex->top_vertices[0][0])/
                               (vertex->left_vertices[1][1] - vertex->top_vertices[0][1]);
                }

            }
        }

        /* condition at first top edge */
        else if (j == vertex->top_vertices[0][1]) {

            /* Just assign nnx1 and we are done */
            nnx1[j] = vertex->top_vertices[0][0];
            break;
        }

        /* spatial marching */
        else {

            nnx1[j] = nnx1[j-1] + delta_x1;
        }

        /* update index */
        j++;

    }

    /* Assign output */
    for (j=0; j<data->ny; j++){
        data->nnx1[j] = floor(nnx1[j]);
    }

    /* clean up */
    free(nnx1);

    return SUCCESS;
}

/* Assemble maximum row-stepping vectors: nnx2 */
int maximum_rowstep(window *data, permuted_vertex *vertex) {

    /* bottom -> right -> top */
    int j=0; 
    double delta_x2;
    double *nnx2;
    const int length = data->ny;
    /* Initialize the nnx2 vector: dimension is ny */
	nnx2 = (double *)malloc(sizeof(double) * length);
    data->nnx2 = (int *)malloc(sizeof(int) * length);
    /* stepping mode */
    data->amode = ROW_STEP;

    /* first make sure the vertices are sorted correctly */
    if (descend_x(vertex->bottom_vertices, vertex->row_size_bv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->bottom_to_right_vertices, vertex->row_size_btrv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->right_vertices, vertex->row_size_rv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->right_to_top_vertices, vertex->row_size_rttv) != SUCCESS) {
            return FAIL;
    }
    if (descend_x(vertex->top_vertices, vertex->row_size_tv) != SUCCESS) {
            return FAIL;
    }

    /* Begin */
    while (j < data->ny) {

        /* Condition at bottom edge */
        if (j == vertex->bottom_vertices[0][1]) {

            /* Assign value to nnx2 */
            nnx2[j] = vertex->bottom_vertices[0][0];

            /* Corner vertex between the bottom edge and right edge */
            if (vertex->bottom_vertices[0][1] == vertex->right_vertices[0][1]) {

                /* No second right vertex */
                if (vertex->row_size_rv < 2) {

                    /* Right-to-top vertices are available */
                    if (vertex->row_size_rttv > 0) {

                        /* Call Hanging_sweep to proceed */
                        j = hanging_sweep(j+1, nnx2, vertex->right_vertices, 0, vertex->top_vertices, 0,
                                          vertex->right_to_top_vertices, vertex->row_size_rttv, data->amode);
                    }

                    /* No right-to-top vertices */
                    else {
                        /* compute delta_x2 using right_vertices[0] and top_vertices[0] */
                        delta_x2 = (vertex->right_vertices[0][0] - vertex->top_vertices[0][0])/
                                   (vertex->right_vertices[0][1] - vertex->top_vertices[0][1]);
                        }

                }
            
                /* Second right vertex exists */
                else {

                    /* Compute delta_x2 using right_vertices[0] and right_vertices[1] */
                    delta_x2 = (vertex->right_vertices[0][0] - vertex->right_vertices[1][0])/
                               (vertex->right_vertices[0][1] - vertex->right_vertices[1][1]);

                }

            
            }

            /* No Corner vertex between the bottom edge and right edge */
            else {
                    
                /* Bottom-to-right vertices exist */
                if (vertex->row_size_btrv > 0) {

                    /* Call hanging sweep to proceed */
                    j = hanging_sweep(j+1, nnx2, vertex->bottom_vertices, 0, vertex->right_vertices, 0,
                                      vertex->bottom_to_right_vertices, vertex->row_size_btrv, data->amode);
                }

                /* Bottom-to-right vertices do not exist */
                else {

                    /* Compute delta_x2 using bottom_vertices[0] and right_vertices[0] */
                    delta_x2 = (vertex->bottom_vertices[0][0] - vertex->right_vertices[0][0])/
                               (vertex->bottom_vertices[0][1] - vertex->right_vertices[0][1]);
                }
            
            }
        }

        /* Condition at first vertex of right edge */
        else if (j == vertex->right_vertices[0][1]) {

            /* Assign value to nnx2 */
            nnx2[j] = vertex->right_vertices[0][0];
            
            /* Second right vertex does not exist */
            if (vertex->row_size_rv < 2) {

                /* Right-to-top vertices exist */
                if (vertex->row_size_rttv > 0) {
                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nnx2, vertex->right_vertices, 0, vertex->top_vertices, 0,
                                      vertex->right_to_top_vertices, vertex->row_size_rttv, data->amode);
                }

                /* Right-to-top vertices do not exist */
                else {
                    /* Compute delta_x2 using right_vertices[0] and top_vertices [0] */
                    delta_x2 = (vertex->right_vertices[0][0] - vertex->top_vertices[0][0])/
                               (vertex->right_vertices[0][1] - vertex->top_vertices[0][1]);
                }

            }
        
            /* Edge propagation: Second right vertex exists */
            else {

                /* Compute delta_x2 using right_vertices[0] and right_vertices[1] */
                delta_x2 = (vertex->right_vertices[0][0] - vertex->right_vertices[1][0])/
                           (vertex->right_vertices[0][1] - vertex->right_vertices[1][1]);
            }
        }

        /* Condition at second vertex of right edge */
        else if (vertex->row_size_rv == 2 && j == vertex->right_vertices[1][1]) {

            /* Corner vertex between the right edge and top edge */
            if (vertex->right_vertices[1][1] == vertex->top_vertices[0][1]) {

                /* Assign value to nnx2 and done */
                nnx2[j] = vertex->right_vertices[1][0];
                break;

            }

            /* No corner vertex betweeen the right edge and top edge*/
            else {
                
                /* Assign value to nnx2 */
                nnx2[j] = vertex->right_vertices[1][0];

                /* Right-to-top hanging vertices exist */
                if (vertex->row_size_rttv > 0) {

                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nnx2, vertex->right_vertices, 1, vertex->top_vertices, 0,
                                      vertex->right_to_top_vertices, vertex->row_size_rttv, data->amode);
                }

                /* No right-to-top vertices */
                else {

                    /* Compute delta_x2 using right_vertics[1] and top_vertices[0] */
                    delta_x2 = (vertex->right_vertices[1][0] - vertex->top_vertices[0][0])/
                               (vertex->right_vertices[1][1] - vertex->top_vertices[0][1]);
                }

            }
        }

        /* Condition at top edge [0] */
        else if (j == vertex->top_vertices[0][1]) {
            /* Assign value */
            nnx2[j] = vertex->top_vertices[0][0];
            break;
        }

        /* spatial marching */
        else {
            nnx2[j] = nnx2[j-1] + delta_x2;
        }

        /* update index */
        j++;
    }

    /* Assign output */
    for (j=0; j<data->ny; j++){
        data->nnx2[j] = ceil(nnx2[j]);
    }

    /* clean up */
    free(nnx2);

    return SUCCESS;
}

/* Assemble minimum column-stepping vectors: nny1 */
int minimum_columnstep(window *data, permuted_vertex *vertex) {

    /* left -> bottom -> right */
    int j=0; 
    double delta_y1;
    double *nny1;
    const int length = data->nx;
    /* Initialize the nny1 vector: dimension is nx */
	nny1 = (double *)malloc(sizeof(double) * length);
    data->nny1 = (int *)malloc(sizeof(int) * length);
    /* stepping mode */
    data->amode = COLUMN_STEP;

    /* first make sure the vertices are sorted correctly */
    if (ascend_y(vertex->left_vertices, vertex->row_size_lv) != SUCCESS) {
            return FAIL;
    }
    if (descend_y(vertex->bottom_to_left_vertices, vertex->row_size_btlv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_x(vertex->bottom_vertices, vertex->row_size_bv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->bottom_to_right_vertices, vertex->row_size_btrv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_y(vertex->right_vertices, vertex->row_size_rv) != SUCCESS) {
            return FAIL;
    }

    while (j < data->nx) {

        /* Condition at first left edge */
        if (j == vertex->left_vertices[0][0]) {

            /* Assign value to nny1 */
            nny1[j] = vertex->left_vertices[0][1];

            /* Corner vertex between the left edge and bottom edge */
            if (vertex->left_vertices[0][0] == vertex->bottom_vertices[0][0]) {

                /* No second bottom vertex */
                if (vertex->row_size_bv < 2) {

                    /* Bottom-to-right vertices exist */
                    if (vertex->row_size_btrv > 0) {

                        /* Call hanging_sweep to proceed */
                        j = hanging_sweep(j+1, nny1, vertex->bottom_vertices, 0, vertex->right_vertices, 0,
                                          vertex->bottom_to_right_vertices, vertex->row_size_btrv, data->amode);
                    }

                    /* No bottom-to-right vertices */
                    else {
                        /* Compute delta_y1 using bottom_vertices[0] and right_vertices[0]*/
                        delta_y1 = (vertex->bottom_vertices[0][1] - vertex->right_vertices[0][1])/
                                   (vertex->bottom_vertices[0][0] - vertex->right_vertices[0][0]);
                    }

                }

                /* Edge propagation: Second bottom vertex exists */
                else {
                    /* Compute delta_y1 using bottom_vertices[0] and bottom_vertices[1] */
                    delta_y1 = (vertex->bottom_vertices[0][1] - vertex->bottom_vertices[1][1])/
                               (vertex->bottom_vertices[0][0] - vertex->bottom_vertices[1][0]);
                }

            }
            /* No corner vertex between the left edge and bottom edge */
            else {

                /* Bottom_to_left hanging vertices exist */
                if (vertex->row_size_btlv > 0) {
                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nny1, vertex->left_vertices, 0, vertex->bottom_vertices, 0,
                                     vertex->bottom_to_left_vertices, vertex->row_size_btlv, data->amode);
                }

                /* No bottom-to-left hanging vertices */
                else {
                    /* Compute delta_y1 using left_vertices[0] and bottom_vertices[0] */
                    delta_y1 = (vertex->left_vertices[0][1] - vertex->bottom_vertices[0][1])/
                               (vertex->left_vertices[0][0] - vertex->bottom_vertices[0][0]);
                }
            }

        }

        /* Condition at first bottom edge */
        else if (j == vertex->bottom_vertices[0][0]) {

            /* Assign value to nny1 */
            nny1[j] = vertex->bottom_vertices[0][1];

            /* No second bottom vertex exists */
            if (vertex->row_size_bv < 2) {

                /* Bottom-to-right vertices exist */
                if (vertex->row_size_btrv > 0) {

                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nny1, vertex->bottom_vertices, 0, vertex->right_vertices, 0,
                                     vertex->bottom_to_right_vertices, vertex->row_size_btrv, data->amode);

                }

                /* No bottom-to-right vertices */
                else {
                    /* Compute delta_y1 using bottom_vertices[0] and right_vertices[0] */
                    delta_y1 = (vertex->bottom_vertices[0][1] - vertex->right_vertices[0][1])/
                               (vertex->bottom_vertices[0][0] - vertex->right_vertices[0][0]);
                }

            }

            /* Edge Propagation: Second bottom verteex exists*/
            else {
                /* Compute delta_y1 using bottom_vertices[0] and bottom_vertices[1] */
                delta_y1 = (vertex->bottom_vertices[0][1] - vertex->bottom_vertices[1][1])/
                           (vertex->bottom_vertices[0][0] - vertex->bottom_vertices[1][0]);
            }

        }

        /* Condition at second bottom edge */
        else if (vertex->row_size_bv == 2 && j == vertex->bottom_vertices[1][0]) {

            /* Corner vertex between the bottom vertices and right vertices */
            if (vertex->bottom_vertices[1][0] == vertex->right_vertices[0][0]) {
                /* Assign value to nny1 and done */
                nny1[j] = vertex->bottom_vertices[1][1];
                break;
            }

            /* No corner vertices */
            else {

                /* Assign value to nny1 */
                nny1[j] = vertex->bottom_vertices[1][1];

                /* Bottom-to-right vertices exists */
                if (vertex->row_size_btrv > 0) {
                    /* Call hanging_vertices to proceed */
                    j = hanging_sweep(j+1, nny1, vertex->bottom_vertices, 1, vertex->right_vertices, 0,
                                     vertex->bottom_to_right_vertices, vertex->row_size_btrv, data->amode);
                }

                /* No bottom-to-right vertices */
                else {
                    /* Compute delta_y1 using bottom_vertices[1] and right_vertices[0] */
                    delta_y1 = (vertex->bottom_vertices[1][1] - vertex->right_vertices[0][1])/
                               (vertex->bottom_vertices[1][0] - vertex->right_vertices[0][0]);
                }

            }

        }

        /* Condition at first right edge: break here */
        else if (j == vertex->right_vertices[0][0]) {

            /* Assign value to nny1 */
            nny1[j] = vertex->right_vertices[0][1];
            break;

        }

        /* Spatial Marching */
        else {
            nny1[j] = nny1[j-1] + delta_y1;
        }

        /* update index */
        j++;
    }

    /* Assign output */
    for (j=0; j<data->nx; j++) {
        data->nny1[j] = floor(nny1[j]);
    }

    /* clean up */
    free(nny1);

    return SUCCESS;
}

/* Assemble maximum column-stepping vectors: nny2 */
int maximum_columnstep(window *data, permuted_vertex *vertex) {

    /* left -> top -> right */
    int j=0; 
    double delta_y2;
    double *nny2;
    const int length = data->nx;
    /* Initialize the nny1 vector: dimension is nx */
	nny2 = (double *)malloc(sizeof(double) * length);
    data->nny2 = (int *)malloc(sizeof(int) * length);
    /* stepping mode */
    data->amode = COLUMN_STEP;

    /* first make sure the vertices are sorted correctly */
    if (descend_y(vertex->left_vertices, vertex->row_size_lv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_x(vertex->left_to_top_vertices, vertex->row_size_lttv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_x(vertex->top_vertices, vertex->row_size_tv) != SUCCESS) {
            return FAIL;
    }
    if (ascend_x(vertex->right_to_top_vertices, vertex->row_size_rttv) != SUCCESS) {
            return FAIL;
    }
    if (descend_y(vertex->right_vertices, vertex->row_size_rv) != SUCCESS) {
            return FAIL;
    }

    while (j < data->nx) { 

        /* Condition at left edge: Already sorted to be max y */
        if (j == vertex->left_vertices[0][0]) {

            /* Assign value */
            nny2[j] = vertex->left_vertices[0][1];

            /* Corner vertex exists between the left edge and top edge */
            if (vertex->left_vertices[0][0] == vertex->top_vertices[0][0]) {

                /* No second top vertex exists */
                if (vertex->row_size_tv < 2) {

                    /* Right_to_top vertices exist */
                    if (vertex->row_size_rttv > 0) {
                        /* Call hanging_sweep to proceed */
                        j = hanging_sweep(j+1, nny2, vertex->top_vertices, 0, vertex->right_vertices, 0,
                                         vertex->right_to_top_vertices, vertex->row_size_rttv, data->amode);

                    }
                    /* No right-to-top vertices */
                    else {
                        /* Compute delta_y2 using top_vertices[0] and right_vertices[0] */
                        delta_y2 = (vertex->top_vertices[0][1] - vertex->right_vertices[0][1])/
                                   (vertex->top_vertices[0][0] - vertex->right_vertices[0][0]);
                    }

                }
                /* Second top vertex exists */
                else {

                    /* Edge propagation: compute delta_y2 using top_vertices[0] and top_vertices[1] */
                    delta_y2 = (vertex->top_vertices[0][1] - vertex->top_vertices[1][1])/
                               (vertex->top_vertices[0][0] - vertex->top_vertices[1][0]);

                }
            }

            /* No corner vertex exists between the left edge and top edge */
            else {

                /* Left-to-top hanging vertices exist */
                if (vertex->row_size_lttv > 0) {
                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nny2, vertex->left_vertices, 0, vertex->top_vertices, 0,
                                      vertex->left_to_top_vertices, vertex->row_size_lttv, data->amode);

                }

                /* No left-to-top hanging vertices */
                else {
                    /* Compute delta_y2 using left_vertices[0] and top_vertices[0]*/
                    delta_y2 = (vertex->left_vertices[0][1] - vertex->top_vertices[0][1])/
                               (vertex->left_vertices[0][0] - vertex->top_vertices[0][0]);
                }

            }

        }

        /* Condition at the first top vertex */
        else if (j == vertex->top_vertices[0][0]) {

            /* Assign value to nny2 */
            nny2[j] = vertex->top_vertices[0][1];

            /* No second top vertices */
            if (vertex->row_size_tv < 2) {

                /* Right-to-top vertices exist */
                if (vertex->row_size_rttv > 0) {
                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nny2, vertex->top_vertices, 0, vertex->right_vertices, 0,
                                      vertex->right_to_top_vertices, vertex->row_size_rttv, data->amode);

                }

                /* No right-to_top vertices */
                else {

                    /* Compute delta_y2 using top_vertices[0] and right_vertices[0] */
                     delta_y2 = (vertex->top_vertices[0][1] - vertex->right_vertices[0][1])/
                                (vertex->top_vertices[0][0] - vertex->right_vertices[0][0]);
                }

            }

            /* Second top vertices exist */
            else {
                /* Edge propagation: Compute delta_y2 using top_vertices[0] and top_vertices[1] */
                delta_y2 = (vertex->top_vertices[0][1] - vertex->top_vertices[1][1])/
                           (vertex->top_vertices[0][0] - vertex->top_vertices[1][0]);

            }

        }

        /* Condition at the second top vertex */
        else if (vertex->row_size_tv == 2 && j == vertex->top_vertices[1][0]) {

            /* Corner vertex between the top edge and right edge */
            if (vertex->top_vertices[1][0] == vertex->right_vertices[0][0] ) {

                /* Assign value and done */
                nny2[j] = vertex->top_vertices[1][1];
                break;

            }

            /* No corner vertex between the top edge and right edge */
            else {

                /* Assign value */
                nny2[j] = vertex->top_vertices[1][1];

                /* Right-to-top vertices exist */
                if (vertex->row_size_rttv > 0) {
                    /* Call hanging_sweep to proceed */
                    j = hanging_sweep(j+1, nny2, vertex->top_vertices, 1, vertex->right_vertices, 0,
                                      vertex->right_to_top_vertices, vertex->row_size_rttv, data->amode);
                }

                /* No right-to-top vertices */
                else {

                    /* Compute delta_y2 using top_vertices[1] and right_vertices[0] */
                    delta_y2 = (vertex->top_vertices[1][1] - vertex->right_vertices[0][1])/
                               (vertex->top_vertices[1][0] - vertex->right_vertices[0][0]);
                }

            }

        }

        /* Condition at right edge: Already sorted to be max y */
        else if (j == vertex->right_vertices[0][0]) {

            /* Assign value */
            nny2[j] = vertex->right_vertices[0][1];
            break;
        }

        /* Spaatial marching */
        else {
            nny2[j] = nny2[j-1] + delta_y2;
        }

        /* update index */
        j++;

    }

    /* Assign output */
    for (j=0; j<data->nx; j++) {
        data->nny2[j] = ceil(nny2[j]);
    }

    /* clean up */
    free(nny2);

    return SUCCESS;
}

/* Assemble the boundary */
int boundary_assembly(window *data, permuted_vertex *vertex) {

    int i;
    const int length_x = data->ny;
    const int length_y = data->nx;
    double **ordered_vertices = NULL;
    int *nnx_diff;
    int *nny_diff;
    nnx_diff = (int *)malloc(sizeof(int) * length_x);
    nny_diff = (int *)malloc(sizeof(int) * length_y);
    if (nnx_diff == NULL || nny_diff == NULL) {
        free(nnx_diff);
        free(nny_diff);
        BLEND_Report(BLEND_MSG_ERROR, "boundary: could not allocate stepping-vector differences\n");
        return FAIL;
    }

    /* Remove duplicate vertices and update row dimensions. TODO: Use pointers here */
    data->row_size = unique_vertices(data->vertices, data->row_size);
    if (data->row_size < 3) {
        free(nnx_diff);
        free(nny_diff);
        return FAIL;
    }

    ordered_vertices = blend_boundary_copy_points(data->vertices, data->row_size);
    if (ordered_vertices == NULL) {
        free(nnx_diff);
        free(nny_diff);
        BLEND_Report(BLEND_MSG_ERROR, "boundary: could not preserve polygon traversal order\n");
        return FAIL;
    }

    /* check that the vertices are bounded */
    if (bounding_boxcheck(data) != SUCCESS) {
        blend_boundary_free_points(ordered_vertices, data->row_size);
        free(nnx_diff);
        free(nny_diff);
        return FAIL;
    }
    blend_boundary_restore_points(data->vertices, ordered_vertices, data->row_size);
    blend_boundary_free_points(ordered_vertices, data->row_size);

    /* Assemble bottom vertices */
    if (assemble_bv(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble left vertices */
    if (assemble_lv(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble top vertices */
    if (assemble_tv(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble right vertices */
    if (assemble_rv(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble traversal-aware interior vertex paths */
    if (assemble_path_vertices(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble minimum row-stepping vectors: nnx1 */
    if (minimum_rowstep(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble maximum row-stepping vectors: nnx2 */
    if (maximum_rowstep(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble minimum column-stepping vectors: nny1 */
    if (minimum_columnstep(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble maximum column-stepping vectors: nny2 */
    if (maximum_columnstep(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Define minnx and minny */
    for (i=0; i<data->ny; i++) {
        nnx_diff[i] = data->nnx2[i] - data->nnx1[i];
    }

    data->minnx = nnx_diff[0];
    for (i=1; i<data->ny; i++){
        if (data->minnx > nnx_diff[i]) {
            data->minnx = nnx_diff[i];
        }
    }
    /* Minimum length of row stepping vectors */
    data->minnx = data->minnx + 1;

    for (i=0; i<data->nx; i++) {
        nny_diff[i] = data->nny2[i] - data->nny1[i];
    }

    data->minny = nny_diff[0];
    for (i=1; i<data->nx; i++){
        if (data->minny > nny_diff[i]) {
            data->minny = nny_diff[i];
        }
    }

    /* Minimum length of column stepping vectors */
    data->minny = data->minny + 1;

    free(nnx_diff);
    free(nny_diff);

    return SUCCESS;
}
