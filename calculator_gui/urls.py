from django.urls import path

from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path("section_a/1", views.kepler_eqn, name="kepler_eqn"),
    path("section_a/2", views.kepler_hyp, name="kepler_hyp"),
    path("section_a/3", views.lagrange, name="lagrange"),
    path("section_b/1", views.section_b_1, name="section b 1"),
    path("section_b/2", views.section_b_2, name="section b 2"),
    path("section_b/3", views.section_b_3, name="section b 3"),
    path("section_b/5", views.section_b_5, name="lambert"),
    path("section_c/1", views.section_c_1, name="inertial to rotating frame"),
    path("section_c/2", views.section_c_2, name="rotating to inertial frame"),
    path("section_c/3", views.section_c_3, name="geocentric to perifocal"),
    path("section_c/4", views.section_c_4, name="perifocal to geocentric"),
    path("section_d/1", views.julian, name="Julian"),
    path("section_d/2", views.lst, name="LST"),
    path("section_f", views.section_f, name="SectionF"),
    path("section_g/1", views.section_g_1, name="Hoffman"),
]
