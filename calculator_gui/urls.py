from django.urls import path

from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path("section_b/1", views.section_b_1, name="section b 1"),
    path("section_b/2", views.section_b_2, name="section b 1"),
    path("section_b/3", views.section_b_3, name="section b 1"),
    path("hoffman", views.hoffman, name="Hoffman"),
    path("julian", views.julian, name="Julian"),
    path("lst", views.lst, name="LST"),
    path("section_a/1", views.kepler_eqn, name="kepler_eqn"),
    path("section_a/2", views.kepler_hyp, name="kepler_hyp"),
    path("section_c/1", views.section_c_1, name="inertial to rotating frame"),
    path("section_c/2", views.section_c_2, name="rotating to inertial frame"),
]
