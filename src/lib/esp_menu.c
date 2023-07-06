/*
 Functions to create, display and manage an ESP-r type menu and abcd...g box
   espmenuinit_ () initialise menu - using list for entries
   updmenu_() initialise popup menu text array
   espmenuitems_ () adds an item to the menu.
   espmenu_ () display menu and wait for selection to be made
   abcbox_callback() actions when an abc... choice is made.
   espabcbox_ () initialise abcbox - will use radio buttons for entries

 */
static const char rcsid[] = "$Id$";

#include <stdio.h>
#include <gtk/gtk.h>
#include <esp_list.h>
#include <gtk/gtkdialog.h>
#include <string.h>
#include "esp-r.h"

/* required by e_menu list and callbacks */
static const gchar *list_item_data_key = "esp_list_item_data";

GtkWidget *gtkplist;	/* the list used for popup menus */
GtkWidget *menu_items;
gint menu_pick;
gint pmenu_pick;

gint abc_pick;
static char choice_list[10][42];        /* character arrays for abcdef boxes */
static char choice_type_list[10]; 
static int choice_list_w[10]; 
static int choice_width = 0;    /* current max choice character width */
static int choice_boxes = 0;    /* current number of abcdef boxes */


/* *************** ESRU choice box text update. *************** */
/*
 This function takes an array of strings from f77 and stores
 it in the static array choice_list for subsequent use by other functions.
*/
/* char choice_list[10][42]; character arrays for choices */
/* char choice_type_list[10]; character array representing choice_list array use */
/* int choice_width = 0; current menu max line length */
/* int choice_boxes = 0; current number of active boxes */
void upd_box_choices_(items,itypes,nitmsptr,iw,len_items)
  char      *items;         /* f77 array of box text strings    */
  char      *itypes;        /* f77 character array (nitmsptr wide)    */
  long int  *nitmsptr;      /* number of choice boxes   */
  long int  *iw;            /* actual max char width in choices    */
  int  len_items;           /* length of string from f77    */
{
  int	i, j, k;
  int	m_line = *nitmsptr;
  int l_m1;

  choice_width = *iw;	/* remember width of choice text */
  choice_boxes = m_line;	/* remember number of choices */
  if(m_line == 0)return;	/* don't bother if no lines */
  strncpy(choice_type_list,itypes,(unsigned int)m_line);	/* copy to static array */

  for(i = 0; i < 10; i++) {	/* clear prior widths...  */
    choice_list_w[i] = 0;
  }
  k = 0;
  for(i = 0; i < choice_boxes; i++) {	/* for each choice...  */
    for(j = 0; j < len_items; j++) {	/* for each character...  */
      choice_list[i][j] = items[k];
      k = k +1;   /* increment for next char in items (a fortran string array does not have
                     nulls between strings in array, it just looks like one long string) */
    }
    choice_list[i][len_items] = '\0';	        /* write terminator          */
    f_to_c_l(choice_list[i],&len_items,&l_m1);  /* find actual length l_m1)  */
    choice_list_w[i] = l_m1;                /* save to choice_list_w */
//    fprintf(stderr,"choice_list %s %d %d %d %d %d %d\n",choice_list[i],i,k,choice_list_w[i],choice_width,len_items,l_m1); 
  }
  return;
}

// initial clear of box choice lines list
//  for ( i = 0; i < 10; i++ ) {
//    strncpy(choice_list[i],
//    "                                         ",41);
//  }

/* **** g_get_esp_item_from_list_cb is the call-back for "selection_changed"
 * signal of the menu selection process */
void g_get_esp_item_from_list_cb ( GtkWidget *a_list,
				   gpointer  func_data )
{
  GList *dlist;
  gchar *buf;

  /* Fetch the selected items of the List, treat as read-only! */
  dlist = GTK_LIST (a_list)->selection; /* listbox widget */

  if (!dlist) {
    esp_selected = buf = NULL; /* If deselected, nullify the strings */
    return;
  }

/* Get the item from the doubly linked list,
 * query the data associated with list_item_data_key.
 * then write it to global (esp_selected) */
  buf = g_object_get_data (G_OBJECT (dlist->data), list_item_data_key);

  got_item_n = gtk_list_child_position ((GtkList*)a_list, dlist->data);
  esp_selected = (gchar*) strdup(buf);
/* debug g_print("you clicked item %d: %s -\n", got_item_n, buf); */
  g_main_loop_quit (menu_loop);
}


/* *** espmenuinit_ () initialise menu - using gtk_list_new for entries *** */
void espmenuinit_ (char *title, int len) {
   gchar *title2,*str,*title3;
   GtkWidget *label, *a_list_item;
   GtkWidget *title_frame;
   gint f_height;	/* pixel height of default font */

   /* DEBUG...
      fprintf(stderr,"TITLE %s\n",title);
      fprintf(stderr,"LEN %d\n",len);
   */

/* create the outer frame 5 pixels wide */
   menu_frame = gtk_vbox_new (FALSE,5);
   gtk_container_set_border_width (GTK_CONTAINER (menu_frame), 1);
   gtk_container_add (GTK_CONTAINER (emenu), menu_frame);
   gtk_widget_show ( menu_frame);

/* create a scolling window for the menu.
 * << still to be done is to set it so that scroll does not
 * << show when there is no need to scroll
 */
   e_scrolled_window = gtk_scrolled_window_new (NULL, NULL);
   gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (e_scrolled_window),
     GTK_POLICY_AUTOMATIC,GTK_POLICY_AUTOMATIC);
   gtk_container_add (GTK_CONTAINER (menu_frame), e_scrolled_window);
   gtk_widget_set_size_request (e_scrolled_window, -1, 500);
   gtk_widget_show (e_scrolled_window);

/* setup the list structure within e_scrolled_window */
   gtklist = gtk_list_new ();
   gtk_scrolled_window_add_with_viewport (GTK_SCROLLED_WINDOW (e_scrolled_window),gtklist);
   gtk_widget_show (gtklist);

   g_signal_connect (G_OBJECT (gtklist), "selection_changed",
                     G_CALLBACK (g_get_esp_item_from_list_cb),NULL);

/* Include the menu title as the first item in the menu. Use g_strndup to
 * go from fortran and then use g_strjoin to prepend a space. */
   title2 = g_strndup(title, (gsize) len);
/*   title3 = g_strjoin(" ","<bold>",title2,"</bold>",NULL);  for when menomonics avail */
   title3 = g_strjoin(" ",title2,NULL);
   label = gtk_label_new (title3);
/*   gtk_label_set_markup_with_mnemonic (GTK_LABEL (label), title3); */
   a_list_item = gtk_list_item_new_with_label (title3);

   gtk_container_add (GTK_CONTAINER (gtklist), a_list_item);
   gtk_widget_show (a_list_item);
   gtk_label_get (GTK_LABEL (label), &str);
   g_object_set_data (G_OBJECT (a_list_item), list_item_data_key, str);
   g_free (title2);

   gtk_widget_show_all (emenu); /* Display the new widgets. */
}


/* *************** ESRU menu text update. *************** */
/*
 This function takes an array of strings from f77 and stores
 it in the static array m_list for subsequent use by other functions.
*/
/* char m_list[MENU_LIST_LEN][125]; character arrays for menu buffer */
/* char mtype_list[MENU_LIST_LEN]; character array representing m_list array use */
/* int m_width = 0; current menu max line length */
/* int m_lines = 0; current number of active menu lines */
void updmenu_(items,itypes,nitmsptr,iw,len_items)
  char      *items;         /* f77 array of menu text strings    */
  char      *itypes;        /* f77 character array (nitmsptr wide)    */
  long int  *nitmsptr;      /* number of menu lines to display   */
  long int  *iw;            /* actual max char width in items    */
  int  len_items;           /* length of menu string from f77    */
{
  int	i, j, k;
  int	pm_line = (int) *nitmsptr;

  pm_width = (int) *iw;	/* remember width of menu text */
  pm_lines = pm_line;	/* remember number of lines */
  if(pm_lines == 0)return;	/* don't bother if no lines */
  strncpy(pmtype_list,itypes,(unsigned int)pm_lines);	/* copy to static array */

  k = 0;
  for(i = 0; i < pm_lines; i++) {	/* for each line...  */
    for(j = 0; j < len_items; j++) {	/* for each character...  */
      pm_list[i][j] = items[k];
      k = k +1;   /* increment for next char in items (a fortran string array does not have
                     nulls between strings in array, it just looks like one long string) */
    }
    pm_list[i][len_items] = '\0';	/* write terminator  */
/* debug  fprintf(stderr,"pm_list %s %d %d %d %d\n",pm_list[i],i,k,pm_width,len_items);  */
  }
  return;
}

/* *** espmenuitems_ () adds an item to the e_menu *** */
void espmenuitems_ (char *item,long int *ino, int len)
{

   GtkWidget *label, *a_list_item;
   gchar *item2,*str;

   item2 = g_strndup(item, (gsize) len);

   /* DEBUG...
      fprintf(stderr,"espmenuitems INO LEN %ld %d %s\n",*ino,len,item2);
    */
   label = gtk_label_new (item2);
   a_list_item = gtk_list_item_new_with_label (item2);

   gtk_container_add (GTK_CONTAINER (gtklist), a_list_item);
   gtk_widget_show (a_list_item);
   gtk_label_get (GTK_LABEL (label), &str);
   g_object_set_data (G_OBJECT (a_list_item), list_item_data_key, str);
   g_free (item2);
   return;
}

/* *** espmenu_ () display menu and wait for selection to be made *** */
void espmenu_ (int *ino)
{

/* Show menu and reset selection. */
   gtk_widget_show_all (gtklist);
   menu_pick = 0;

   /* DEBUG...
    * fprintf(stderr,"INO1 %d\n",*ino);
    * fprintf(stderr,"INO2 %d\n",menu_pick);
    */

/* Create a new event loop and run it.  When a selection is
 * made the value of got_item_n will equal to the position in
 * the list (because the menu title is the first item).
 * This value is placed to the location that ino points to.
 */
   g_main_loop_run (menu_loop);
   *ino = got_item_n;	/* got_item_n is global set by g_get_esp_item_from_list_cb */

   /* DEBUG...
    * fprintf(stderr,"INO4 %d\n",*ino);
    * fprintf(stderr,"INO5 %d\n",got_item_n);
    */
   /* get rid of the entities used in the previous dynamic menu */
    gtk_widget_destroy (gtklist);
    gtk_widget_destroy (e_scrolled_window);
    gtk_widget_destroy (menu_frame);
}


/* *** abcbox_callback() actions when an abc... choice is made. *** */

void abcbox_callback( GtkWidget *widget,
                      gpointer  data )
{
/* debug g_print ("Hello again - %d was pressed\n", GPOINTER_TO_INT (data)); */
  abc_pick = GPOINTER_TO_INT (data);	/* see spinbutton example */
}

/* *** espabcbox_ () initialise abcbox - will use radio buttons for entries *** */
void espabcbox_ (char *msg1, char *aopt, char *bopt, char *copt,
                 char *dopt, char *eopt, char *fopt, char *gopt, long int *ipick,
                 int msg1_len, int aopt_len, int bopt_len, int copt_len,
                 int dopt_len, int eopt_len, int fopt_len, int gopt_len) {

   GtkWidget *askbox, *hbox, *left_col, *right_col, *button, *label;
   GSList *group;
   gchar *title_local, *msg1_local;
   gchar *aopt_local, *bopt_local, *copt_local, *dopt_local, *eopt_local, *fopt_local;
   gchar *gopt_local;
   gchar *question_local;
   gint result;
   gint aopt_l,bopt_l,copt_l,dopt_l,eopt_l,fopt_l,gopt_l; /* non-blank lengths for options */
   int msg1_l; /* non-blank lengths for prompt */

   title_local = "  ";

/* find out actual length of each prompt and then total length with a space between. */
   f_to_c_l(msg1,&msg1_len,&msg1_l);
   question_local = g_strndup(msg1, (gsize) msg1_l);
/* debug  fprintf(stderr,"ask phrase %s\n",question_local); */

   f_to_c_l(aopt,&aopt_len,&aopt_l);
   aopt_local = g_strndup(aopt, (gsize) aopt_l);
   f_to_c_l(bopt,&bopt_len,&bopt_l);
   bopt_local = g_strndup(bopt, (gsize) bopt_l);
   f_to_c_l(copt,&copt_len,&copt_l);
   copt_local = g_strndup(copt, (gsize) copt_l);
   f_to_c_l(dopt,&dopt_len,&dopt_l);
   dopt_local = g_strndup(dopt, (gsize) dopt_l);
   f_to_c_l(eopt,&eopt_len,&eopt_l);
   eopt_local = g_strndup(eopt, (gsize) eopt_l);
   f_to_c_l(fopt,&fopt_len,&fopt_l);
   fopt_local = g_strndup(fopt, (gsize) fopt_l);
   f_to_c_l(gopt,&gopt_len,&gopt_l);
   gopt_local = g_strndup(gopt, (gsize) gopt_l);
/* debug fprintf(stderr,"non-blank lengths are %d %d %d %d %d %d %d\n",aopt_l,bopt_l,copt_l,dopt_l,eopt_l,fopt_l,gopt_l); */

   /* Set ok response, but if *ipick is zero reset to one. */
   abc_pick = (gint) *ipick;
   if ( abc_pick == 0 ) abc_pick = 1;

   /* Create the widgets: first the dialog window, then split the default
   vbox into two (left_col and right_col) by using a hbox so as the items
   can be displayed in two columns.
   */
   askbox = gtk_dialog_new_with_buttons(title_local,
     GTK_WINDOW (window),GTK_DIALOG_DESTROY_WITH_PARENT,
     GTK_STOCK_HELP, GTK_RESPONSE_HELP,GTK_STOCK_OK, GTK_RESPONSE_OK,
     GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,NULL);
   label = gtk_label_new (question_local);
   gtk_container_add (GTK_CONTAINER (GTK_DIALOG(askbox)->vbox), label);

   hbox = gtk_hbox_new (TRUE, 2);
   gtk_container_add (GTK_CONTAINER (GTK_DIALOG(askbox)->vbox),hbox);
   left_col = gtk_vbox_new (TRUE, 2);
   gtk_box_pack_start (GTK_BOX (hbox),left_col, TRUE, TRUE, 0);
   right_col = gtk_vbox_new (TRUE, 2);
   gtk_box_pack_start (GTK_BOX (hbox),right_col, TRUE, TRUE, 0);

   /* Add entries to the columns in turn */
    button = gtk_radio_button_new_with_label (NULL, aopt_local);
    gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
    if ( abc_pick==1 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
    g_signal_connect (G_OBJECT (button), "pressed",
                      G_CALLBACK (abcbox_callback), GINT_TO_POINTER (1));
    gtk_widget_show (button);
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));

    button = gtk_radio_button_new_with_label (group, bopt_local);
    gtk_box_pack_start (GTK_BOX (right_col), button, TRUE, TRUE, 0);
    if ( abc_pick==2 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
    g_signal_connect (G_OBJECT (button), "pressed",
                      G_CALLBACK (abcbox_callback), GINT_TO_POINTER (2));
    gtk_widget_show (button);

    /* Additional entries.  Check if defined before displaying (string longer than 1 char) */
    if (copt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), copt_local);
      gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
      if ( abc_pick==3 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (3));
      gtk_widget_show (button);
    }

    if (dopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), dopt_local);
      gtk_box_pack_start (GTK_BOX (right_col), button, TRUE, TRUE, 0);
      if ( abc_pick==4 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (4));
      gtk_widget_show (button);
    }

    if (eopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), eopt_local);
      gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
      if ( abc_pick==5 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                      G_CALLBACK (abcbox_callback), GINT_TO_POINTER (5));
      gtk_widget_show (button);
    }

    if (fopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), fopt_local);
      gtk_box_pack_start (GTK_BOX (right_col), button, TRUE, TRUE, 0);
      if ( abc_pick==6 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (6));
      gtk_widget_show (button);
    }

    if (gopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), gopt_local);
      gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
      if ( abc_pick==7 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (7));
      gtk_widget_show (button);
    }

   /*
      Display the new widgets.
   */
   gtk_widget_show_all (askbox);

   /* Set dialog properties and wait for user response */
   gtk_window_set_modal (GTK_WINDOW (askbox), TRUE);
   gtk_window_set_transient_for(GTK_WINDOW (askbox), GTK_WINDOW (window));

   result = gtk_dialog_run (GTK_DIALOG (askbox));
   switch (result)
      {
       case GTK_RESPONSE_OK:
       /*   fprintf(stderr,"Goodbye - %d %d was selected\n", abc_pick,result); debug */
          *ipick = (long int) abc_pick;
          break;
       case GTK_RESPONSE_CANCEL:
       /*   fprintf(stderr,"Goodbye - %d %d was original choice\n", abc_pick,result); debug */
          *ipick = -3;
          break;
       case GTK_RESPONSE_HELP:
       /*   fprintf(stderr,"Goodbye - %d %d help was requested\n", abc_pick,result); debug */
          *ipick = -8;
          break;
       default:
       /*   fprintf(stderr,"Goodbye - %d %d fell to default response\n", abc_pick,result); debug */
          *ipick = -1;
          break;
      }
   gtk_widget_destroy (askbox);

}


/* *** espmbox_ () initialise abcbox - will use radio buttons for entries *** */
void espmbox_ (char *msg1, long int *ipick, int msg1_len) {

   GtkWidget *askbox, *hbox, *left_col, *right_col, *button, *label;
   GSList *group;
   gchar *title_local, *msg1_local;
   gchar *aopt_local, *bopt_local, *copt_local, *dopt_local, *eopt_local, *fopt_local;
   gchar *gopt_local;
   gchar *question_local;
   gint result;
   gint aopt_l,bopt_l,copt_l,dopt_l,eopt_l,fopt_l,gopt_l; /* non-blank lengths for options */
   int msg1_l; /* non-blank lengths for prompt */
   int nopts;	/* number of options (based on if option text blank) */

   title_local = "  ";

/* find out actual length of each prompt and then total length with a space between. */
   f_to_c_l(msg1,&msg1_len,&msg1_l);
   question_local = g_strndup(msg1, (gsize) msg1_l);
/* debug  fprintf(stderr,"ask phrase %s\n",question_local); */

   nopts = 0;
   aopt_l=choice_list_w[0];
   aopt_local = g_strndup(choice_list[0], (gsize) aopt_l);
   if ( choice_list_w[1] > 1 ) nopts = 2;
   bopt_l=choice_list_w[1];
   bopt_local = g_strndup(choice_list[1], (gsize) bopt_l);
   if ( choice_list_w[2] > 1 ) nopts = 3;
   copt_l=choice_list_w[2];
   copt_local = g_strndup(choice_list[2], (gsize) copt_l);
   if ( choice_list_w[3] > 1 ) nopts = 4;
   dopt_l=choice_list_w[3];
   dopt_local = g_strndup(choice_list[3], (gsize) dopt_l);
   if ( choice_list_w[4] > 1 ) nopts = 5;
   eopt_l=choice_list_w[4];
   eopt_local = g_strndup(choice_list[4], (gsize) eopt_l);
   if ( choice_list_w[5] > 1 ) nopts = 6;
   fopt_l=choice_list_w[5];
   fopt_local = g_strndup(choice_list[5], (gsize) fopt_l);
   if ( choice_list_w[6] > 1 ) nopts = 7;
   gopt_l=choice_list_w[6];
   gopt_local = g_strndup(choice_list[6], (gsize) gopt_l);
// fprintf(stderr,"non-blank lengths are %d %d %d %d %d %d %d\n",aopt_l,bopt_l,copt_l,dopt_l,eopt_l,fopt_l,gopt_l);

   /* Set ok response, but if *ipick is zero reset to one. */
   abc_pick = (gint) *ipick;
   if ( abc_pick == 0 ) abc_pick = 1;

   /* Create the widgets: first the dialog window, then split the default
   vbox into two (left_col and right_col) by using a hbox so as the items
   can be displayed in two columns.
   */
   askbox = gtk_dialog_new_with_buttons(title_local,
     GTK_WINDOW (window),GTK_DIALOG_DESTROY_WITH_PARENT,
     GTK_STOCK_HELP, GTK_RESPONSE_HELP,GTK_STOCK_OK, GTK_RESPONSE_OK,
     GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,NULL);
   label = gtk_label_new (question_local);
   gtk_container_add (GTK_CONTAINER (GTK_DIALOG(askbox)->vbox), label);

   hbox = gtk_hbox_new (TRUE, 2);
   gtk_container_add (GTK_CONTAINER (GTK_DIALOG(askbox)->vbox),hbox);
   left_col = gtk_vbox_new (TRUE, 2);
   gtk_box_pack_start (GTK_BOX (hbox),left_col, TRUE, TRUE, 0);
   right_col = gtk_vbox_new (TRUE, 2);
   gtk_box_pack_start (GTK_BOX (hbox),right_col, TRUE, TRUE, 0);

   /* Add entries to the columns in turn */
    button = gtk_radio_button_new_with_label (NULL, aopt_local);
    gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
    if ( abc_pick==1 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
    g_signal_connect (G_OBJECT (button), "pressed",
                      G_CALLBACK (abcbox_callback), GINT_TO_POINTER (1));
    gtk_widget_show (button);
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));

    button = gtk_radio_button_new_with_label (group, bopt_local);
    gtk_box_pack_start (GTK_BOX (right_col), button, TRUE, TRUE, 0);
    if ( abc_pick==2 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
    g_signal_connect (G_OBJECT (button), "pressed",
                      G_CALLBACK (abcbox_callback), GINT_TO_POINTER (2));
    gtk_widget_show (button);

    /* Additional entries.  Check if defined before displaying (string longer than 1 char) */
    if (copt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), copt_local);
      gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
      if ( abc_pick==3 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (3));
      gtk_widget_show (button);
    }

    if (dopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), dopt_local);
      gtk_box_pack_start (GTK_BOX (right_col), button, TRUE, TRUE, 0);
      if ( abc_pick==4 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (4));
      gtk_widget_show (button);
    }

    if (eopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), eopt_local);
      gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
      if ( abc_pick==5 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                      G_CALLBACK (abcbox_callback), GINT_TO_POINTER (5));
      gtk_widget_show (button);
    }

    if (fopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), fopt_local);
      gtk_box_pack_start (GTK_BOX (right_col), button, TRUE, TRUE, 0);
      if ( abc_pick==6 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (6));
      gtk_widget_show (button);
    }

    if (gopt_l > 1 ) {
      button = gtk_radio_button_new_with_label_from_widget
                            (GTK_RADIO_BUTTON (button), gopt_local);
      gtk_box_pack_start (GTK_BOX (left_col), button, TRUE, TRUE, 0);
      if ( abc_pick==7 ) {gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);}
      g_signal_connect (G_OBJECT (button), "pressed",
                        G_CALLBACK (abcbox_callback), GINT_TO_POINTER (7));
      gtk_widget_show (button);
    }

   /*
      Display the new widgets.
   */
   gtk_widget_show_all (askbox);

   /* Set dialog properties and wait for user response */
   gtk_window_set_modal (GTK_WINDOW (askbox), TRUE);
   gtk_window_set_transient_for(GTK_WINDOW (askbox), GTK_WINDOW (window));

   result = gtk_dialog_run (GTK_DIALOG (askbox));
   switch (result)
      {
       case GTK_RESPONSE_OK:
          /* fprintf(stderr,"Goodbye - %d %d was selected\n", abc_pick,result); debug */
          *ipick = (long int) abc_pick;
          break;
       case GTK_RESPONSE_CANCEL:
          /* fprintf(stderr,"Goodbye - %d %d was original choice\n", abc_pick,result); debug */
          *ipick = -3;
          break;
       case GTK_RESPONSE_HELP:
          /* fprintf(stderr,"Goodbye - %d %d help was requested\n", abc_pick,result); debug */
          *ipick = -8;
          break;
       default:
          /* fprintf(stderr,"Goodbye - %d %d fell to default response\n", abc_pick,result); debug */
          *ipick = -1;
          break;
      }
   gtk_widget_destroy (askbox);

}


