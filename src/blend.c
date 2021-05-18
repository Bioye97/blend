#include "blend.h"

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
        fprintf(stderr, "The number of unique vertices is less than three. Please revise vertices. \n");
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

int bounding_boxcheck(window_t *data) {

    /* Ensure that the range of the vertices is 
    [0 nx-1] and [0 ny-1] */

    /* Check that the polygon represented by the vertices 
    are bounded by the model dimension */

    /* Minimum bound check on x */
    if (ascend_x(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][0] < 0) {
        fprintf(stderr, "Minimum value for the x-grid coordinates is less than 0. \n");
        fprintf(stderr, "A vertex is outside the model domain. Please revise vertices. \n");
        return FAIL;
    }

    if (data->vertices[0][0] > 0) {
        fprintf(stderr, "Minimum value for the x-grid coordinates is greater than 0. \n");
        fprintf(stderr, "The polygon is not bounded by the model domain. Please revise vertices. \n");
        return FAIL;
    }

    /* Maximum bound check on x */
     if (descend_x(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][0] < data->nx - 1) {
        fprintf(stderr, "Maximum value for the x-grid coordinates is less than %d. \n", data->nx-1);
        fprintf(stderr, "The polygon is not bounded by the model domain. Please revise vertices. \n");
        return FAIL;
    }

    if (data->vertices[0][0] > data->nx - 1) {
        fprintf(stderr, "Maximum value for the x-grid coordinates is greater than %d. \n", data->nx-1);
        fprintf(stderr, "A vertex is outside the model domain. Please revise vertices. \n");
        return FAIL;
    }

    /* Minimum bound check on y */
    if (ascend_y(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][1] < 0) {
        fprintf(stderr, "Minimum value for the y-grid coordinates is less than 0. \n");
        fprintf(stderr, "A vertex is outside the model domain. Please revise vertices. \n");
        return FAIL;
    }

    if (data->vertices[0][1] > 0) {
        fprintf(stderr, "Minimum value for the y-grid coordinates is greater than 0. \n");
        fprintf(stderr, "The polygon is not bounded by the model domain. Please revise vertices. \n");
        return FAIL;
    }

    /* Maximum bound check on y */
    if (descend_y(data->vertices, data->row_size) != SUCCESS) {
        return FAIL;
    }

    if (data->vertices[0][1] < data->ny - 1) {
        fprintf(stderr, "Maximum value for the y-grid coordinates is less than %d. \n", data->ny-1);
        fprintf(stderr, "The polygon is not bounded by the model domain. Please revise vertices. \n");
        return FAIL;
    }

    if (data->vertices[0][1] > data->ny - 1) {
        fprintf(stderr, "Maximum value for the y-grid coordinates is greater than %d. \n", data->ny-1);
        fprintf(stderr, "A vertex is outside the model domain. Please revise vertices. \n");
        return FAIL;
    }

    return SUCCESS;
}

/* Assemble bottom vertices */
int assemble_bv(window_t *data, permuted_vertex_t *vertex) {

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
       fprintf(stderr, "No bottom vertices found. Please revise vertices. \n");
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble left vertices */
int assemble_lv(window_t *data, permuted_vertex_t *vertex) {
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
       fprintf(stderr, "No left vertices found. Please revise vertices. \n");
    }

    free(tmp_array);
    return SUCCESS;

}

/* Assemble top vertices */
int assemble_tv(window_t *data, permuted_vertex_t *vertex) {
    
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
       fprintf(stderr, "No top vertices found. Please revise vertices. \n");
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble right vertices */
int assemble_rv(window_t *data, permuted_vertex_t *vertex) {
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
       fprintf(stderr, "No right vertices found. Please revise vertices. \n");
    }

    free(tmp_array);
    return SUCCESS;
}

/* Assemble bottom-to-left-vertices */
int assemble_btlv(window_t *data, permuted_vertex_t *vertex) {

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
                fprintf(stderr, "Irregular vertices detected from the bottom edge \
                to the left edge of the bounding box. Revise vertices in this range \
                as the y coordinates should be strictly increasing. \n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] < tmp_array[m+1][0]) {
                fprintf(stderr, "Irregular vertices detected from the bottom edge \
                to the left edge of the bounding box. Revise vertices in this range \
                as the x coordinates should be strictly decreasing. \n");
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
int assemble_lttv(window_t *data, permuted_vertex_t *vertex){

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
                fprintf(stderr, "Irregular vertices detected from the left edge \
                to the top edge of the bounding box. Revise vertices in this range \
                as the y coordinates should be strictly increasing. \n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] > tmp_array[m+1][0]) {
                fprintf(stderr, "Irregular vertices detected from the left edge \
                to the top edge of the bounding box. Revise vertices in this range \
                as the x coordinates should be strictly increasing. \n");
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
int assemble_rttv(window_t *data, permuted_vertex_t *vertex) {

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
                fprintf(stderr, "Irregular vertices detected from the right edge \
                to the top edge of the bounding box. Revise vertices in this range \
                as the y coordinates should be strictly increasing. \n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] < tmp_array[m+1][0]) {
                fprintf(stderr, "Irregular vertices detected from the right edge \
                to the top edge of the bounding box. Revise vertices in this range \
                as the x coordinates should be strictly decreasing. \n");
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
int assemble_btrv(window_t *data, permuted_vertex_t *vertex) {
        
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
                fprintf(stderr, "Irregular vertices detected from the bottom edge \
                to the right edge of the bounding box. Revise vertices in this range \
                as the y coordinates should be strictly increasing. \n");
                return FAIL;
            }
        }
        for (m=0; m<count-1; m++) {
            if (tmp_array[m][0] > tmp_array[m+1][0]) {
                fprintf(stderr, "Irregular vertices detected from the bottom edge \
                to the right edge of the bounding box. Revise vertices in this range \
                as the x coordinates should be strictly increasing. \n");
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

/* Hanging sweep */
int hanging_sweep(int j, double *vector, double **prev_array, int prev_row, double **sub_array, int sub_row,
                  double **hang_array, int hang_len,  assembly_type_t step_type) {

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
                index--; // update index for next vertex condition
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
                index--; // update index for next vertex condition
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
int minimum_rowstep(window_t *data, permuted_vertex_t *vertex) {

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
int maximum_rowstep(window_t *data, permuted_vertex_t *vertex) {

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
int minimum_columnstep(window_t *data, permuted_vertex_t *vertex) {

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
int maximum_columnstep(window_t *data, permuted_vertex_t *vertex) {

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
int boundary_assembly(window_t *data, permuted_vertex_t *vertex) {

    int i;
    const int length_x = data->ny;
    const int length_y = data->nx;
    int *nnx_diff;
    int *nny_diff;
    nnx_diff = (int *)malloc(sizeof(int) * length_x);
    nny_diff = (int *)malloc(sizeof(int) * length_y);

    /* Remove duplicate vertices and update row dimensions. TODO: Use pointers here */
    data->row_size = unique_vertices(data->vertices, data->row_size);

    /* check that the vertices are bounded */
    if (bounding_boxcheck(data) != SUCCESS) {
        return FAIL;
    }

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

    /* Assemble bottom-to-left-vertices */
    if (assemble_btlv(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble left-to-top vertices */
    if (assemble_lttv(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble right-to-top vertices */
    if (assemble_rttv(data, vertex) != SUCCESS) {
        return FAIL;
    }

    /* Assemble bottom-to-right vertices */
    if (assemble_btrv(data, vertex) != SUCCESS) {
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

/* boxcar function */
double boxcar(void) {
    return 1.0;
}

/* cosine function */
double cosine(int x, int n, int nmax, int nmin, double r1, double r2) {
    double ni, nf;
    double value, t;

    ni = r1 * (double)nmax;
    nf = r2 * (double)nmax;

    if (ni > floor(nmin/2)) {
        ni = r1 * (double)n;
        nf = r2 * (double)n;
    }

    if (x > ni && x <= n - nf ) {
        value = 1;
    }

    else if (x <= ni) {
        t = ((2 * (double)x) - 1)/(2 * (double)ni);
        value = 0.5 * (1 - cos(M_PI * t));
    }

    else {

        t = (2 * ((double)x - (double)n + 2 * (double)nf) - 1)/(2 * (double)nf);
        value = 0.5 * (1 - cos(M_PI * t));
    }

    return value;
}

/* scosine function */
double scosine(int x, int n, int nmax, int nmin, double r) {
    double ni;
    double t, value;

    ni = r * (double)nmax;

    if (ni > floor(nmin/2)) {
        ni = r * (double)n;
    }

    if (x > ni) {
        value = 1;
    }
    else {
        t = ((2 * (double)x) - 1)/(2 * (double)ni);
        value = 0.5 * (1 - cos(M_PI * t));
    }

    return value;

}

/* ecosine function */
double ecosine(int x, int n, int nmax, int nmin, double r) {

    double nf;
    double value, t;

    nf = r * (double)nmax;

    if (nf > floor(nmin/2)) {
        nf = r * (double)n;
    }

    if (x > n - nf) {

        t = (2 * ((double)x - (double)n + 2 * (double)nf) - 1)/(2 * (double)nf);
        value = 0.5 * (1 - cos(M_PI * t));
    }
    else {
        value = 1;
    }

    return value;
}

/* trapezoid function */
double trapezoid(int x, int n, int nmax, int nmin, double r1, double r2) {
    double ni, nf;
    double value, t;
    ni = r1 * (double)nmax;
    nf = r2 * (double)nmax;

    if (ni > floor(nmin/2)) {
        ni = r1 * (double)n;
        nf = r2 * (double)n;
    }

    if (x > ni && x <= n - nf ) {
        value = 1;
    }

    else if (x <= ni) {
        t = ((2 * (double)x) - 1)/(2 * (double)ni);
        value = t;
    }

    else {
        t = (2 * ((double)x - (double)n + 2 * (double)nf) - 1)/(2 * (double)nf);
        value = 2 - t;
    }

    return value;
}

/* strapezoid function */
double strapezoid(int x, int n, int nmax, int nmin, double r) {
    double ni;
    double t, value;

    ni = r * (double)nmax;

    if (ni > floor(nmin/2)) {
        ni = r * (double)n;
    }

    if (x > ni) {
        value = 1;
    }
    else {
        t = ((2 * (double)x) - 1)/(2 * (double)ni);
        value = t;
    }

    return value;

}

/* etrapezoid function */
double etrapezoid(int x, int n, int nmax, int nmin, double r) {

    double nf;
    double value, t;

    nf = r * (double)nmax;

    if (nf > floor(nmin/2)) {
        nf = r * (double)n;
    }

    if (x > n - nf) {

        t = (2 * ((double)x - (double)n + 2 * (double)nf) - 1)/(2 * (double)nf);
        value = 2 - t;
    }
    else {
        value = 1;
    }

    return value;
}

/* Smoothing function */
double window_function(int x, int x1, int x2, int n, int nmin, double r, char func[MAX_WINDOWFUNC_LEN]) {

    int n1, n2;
    double taper;

    if (r >= 0.5) {
        fprintf(stderr, "Taper ratio cannot exceed or be equal to 0.5. \n");
        return WFUNC_ERROR;
    }

    if (x < 1 || x > n) {
        fprintf(stderr, "Query point is outside model domain. \n");
        return WFUNC_ERROR;
    }

    /* get lengths of the first and second segments */
    n1 = x1 - 1;
    n2 = x2 - x1 + 1;

    /* Check if the query point lies outside the polygon */
    if (x < x1){

        return 0;
    }

    else if (x > x2) {
        
        return 0;
    }

    else {

        /* First remap query point to new domain */
        x = x - n1;

        /* return taper using appropriate function */
        if (strcmp(func, WFUNC_BOXCAR) == 0) {
            taper = boxcar();
        }

        else if (strcmp(func, WFUNC_COSINE) == 0) {
            
            taper = cosine(x, n2, n, nmin, r, r);
        }

        else if (strcmp(func, WFUNC_COSINE_START) == 0) {
            taper = scosine(x, n2, n, nmin, r);
        }

        else if (strcmp(func, WFUNC_COSINE_END) == 0) {
            taper = ecosine(x, n2, n, nmin, r);
        }

         else if (strcmp(func, WFUNC_TRAPEZOID) == 0) {
            taper = trapezoid(x, n2, n, nmin, r, r);
        }

         else if (strcmp(func, WFUNC_TRAPEZOID_START) == 0) {
            taper = strapezoid(x, n2, n, nmin, r);
        }

         else if (strcmp(func, WFUNC_TRAPEZOID_END) == 0) {
            taper = etrapezoid(x, n2, n, nmin, r);
        }
        else {
            fprintf(stderr, "The window function you specified does not exist. \n");
            return WFUNC_ERROR;
        }
    }
    return taper;

}

/* Embedding contribution */
int embedding_contribution(int x, int y, int z, window_t *data) {

    int nny1, nny2, nnx1, nnx2;
    double taper_x, taper_y, taper_z;
    /* First define the values of these vectors */
    nny1 = data->nny1[x];
    nny2 = data->nny2[x];
    nnx1 = data->nnx1[y];
    nnx2 = data->nnx2[y];

    /* get taper_x: TODO: remove the +1 after generalization  */
    taper_x = window_function(y+1, nny1+1, nny2+1, data->ny, data->minny, data->ratio, data->x_function);
    /* Error check */
    if (taper_x == WFUNC_ERROR) {
        fprintf(stderr, "An error occured when trying to compute the weight from x_function. \n");
        return FAIL;
    }


    /* get taper_y: TODO: remove the +1 after generalization */
    taper_y = window_function(x+1, nnx1+1, nnx2+1, data->nx, data->minnx, data->ratio, data->y_function);
    /* Error check */
    if (taper_y == WFUNC_ERROR) {
        fprintf(stderr, "An error occured when trying to compute the weight from y_function. \n");
        return FAIL;
    }

    /* get taper_z TODO: remove the +1 after generalization */
    taper_z = window_function(z+1, 1, data->nz, data->nz, data->nz, data->ratio, data->z_function);
    /* Error check */
    if (taper_z == WFUNC_ERROR) {
        fprintf(stderr, "An error occured when trying to compute the weight from z_function. \n");
        return FAIL;
    }

    /* Assemble output */
    data->contribution = taper_x * taper_y * taper_z;

    return SUCCESS;
}